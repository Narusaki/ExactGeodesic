// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the ICH_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// ICH_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef ICH_EXPORTS
#define ICH_API __declspec(dllexport)
#else
#define ICH_API __declspec(dllimport)
#endif

#include <Mesh.h>
#include <vector>
#include <queue>
#include <math.h>

#define RELATIVE_ERROR 1e-2
//#define L_RELATIVE_ERROR 1e-3

// This class is exported from the ICH.dll
class ICH_API ICH {
public:
	struct Window
	{
		unsigned edgeID;
		double b0, b1, d0, d1;
		double pseudoSrcDist, minDist;
		unsigned srcID, pseudoSrcId;
		int pseudoSrcBirthTime;
		int level;

		bool operator< (const Window& right) const
		{
			return minDist > right.minDist;
		}

		void calcMinDist()
		{
			double wLen = b1 - b0;
			double xProj = (d0*d0 + wLen*wLen - d1*d1) / (2 * wLen);
			if (xProj < 0.0) minDist = d0 + pseudoSrcDist;
			else if (xProj > wLen) minDist = d1 + pseudoSrcDist;
			else {
				minDist = sqrt(fabs(d0*d0 - xProj*xProj)) + pseudoSrcDist;
			}
		}

		Vector2D FlatenedSrc() const
		{
			Vector2D src2D;
			double wLen = b1 - b0;
			src2D.x = (d0*d0 + wLen*wLen - d1*d1) / (2.0 * wLen);
			src2D.y = sqrt(fabs(d0*d0 - src2D.x*src2D.x));
			src2D.x += b0;
			return src2D;
		}
	};

	struct PseudoWindow
	{
		unsigned vertID;
		double dist;
		unsigned srcId, pseudoSrcId;
		unsigned pseudoBirthTime;
		unsigned level;

		bool operator< (const PseudoWindow &right) const {
			return dist > right.dist;
		};
	};

	struct SplitInfo
	{
		double dist;
		unsigned pseudoSrcId, srcId;
		unsigned level;
		double x;

		SplitInfo() { dist = DBL_MAX; pseudoSrcId = -1; x = DBL_MAX; }
	};

	struct VertInfo
	{
		int birthTime;
		double dist;
		int enterEdge;

		VertInfo() { birthTime = -1; dist = DBL_MAX; enterEdge = -1; }
	};

	struct GeodesicKeyPoint
	{
		bool isVertex;
		unsigned id;
		double pos;
	};
	
public:
	ICH();
	~ICH();
	// TODO: add your methods here.
	void Clear();
	void AssignMesh(CMesh *mesh_);
	void AddSource(unsigned vertId);
	void AddSource(unsigned faceId, Vector3D pos);
	void AddFacesKeptWindow(unsigned faceId);
	void Execute(int totalCalcVertNum_ = -1);
	void OutputStatisticInfo();
	std::list<GeodesicKeyPoint> BuildGeodesicPathTo(unsigned faceId, Vector3D pos, unsigned &srcId);
	std::list<GeodesicKeyPoint> BuildGeodesicPathTo(unsigned vertId, unsigned &srcId);
	double GetDistanceTo(unsigned vertId);

private:
	void Initialize();
	void PropagateWindow(const Window &win);

	void GenSubWinsForPseudoSrc(const PseudoWindow &pseudoWin);
	void GenSubWinsForPseudoSrcFromWindow(const PseudoWindow &pseudoWin, unsigned &startEdge, unsigned &endEdge);
	void GenSubWinsForPseudoSrcFromPseudoSrc(const PseudoWindow &pseudoWin, unsigned &startEdge, unsigned &endEdge);

	bool IsValidWindow(const Window &win, bool isLeftChild);
	void BuildWindow(const Window &fatherWin, 
		unsigned edge, 
		double t0, double t1, 
		const Vector2D &v0, const Vector2D &v1, 
		Window &win);

	double Intersect(const Vector2D &v0, const Vector2D &v1, const Vector2D &p0, const Vector2D &p1);

private:
	CMesh *mesh;
	std::vector< SplitInfo > splitInfos;
	std::vector< VertInfo > vertInfos;
	std::priority_queue< Window > winQ;
	std::priority_queue< PseudoWindow > pseudoSrcQ;
	std::vector< unsigned > sourceVerts;
	std::vector< std::pair< unsigned, Vector3D > > sourcePoints;

	std::vector< Window > storedWindows;
	std::vector< unsigned > keptFaces;

	// statistics
	int numOfWinGen;
	int maxWinQSize, maxPseudoQSize;
	int totalCalcVertNum;
	
};
