#include <Mesh.h>
#include <ICH.h>
#include <fstream>
#include <ctime>
#include <vector>

using namespace std;

vector< int > srcVerts, dstVerts;
vector< pair<int, Vector3D> > srcPoints, dstPoints;

vector<string> splitString(string s, string sep)
{
	vector<string> res;
	if (s.size() == 0) return res;
	size_t pos = s.find(sep);

	while (pos != string::npos)
	{
		if (s.size() == 0) break;
		res.push_back(s.substr(0, pos));
		s = s.substr(pos + sep.size(), s.length() - pos - sep.size());
		pos = s.find(sep);
	}
	if (s.size() != 0) res.push_back(s);
	return res;
}

bool LoadInputFile(const char* fileName)
{
	string curLine;
	ifstream input(fileName);
	if (!input)
	{
		cout << "Cannot open inputfile " << fileName << "!" << endl;
		return false;
	}
	while (getline(input, curLine))
	{
		auto parts = splitString(curLine, " ");
		if (parts[0] == "srcVert")
		{
			for (int i = 1; i < parts.size(); ++i)
				srcVerts.push_back(atoi(parts[i].c_str()));
		}
		else if (parts[0] == "srcPoint")
		{
			for (int i = 1; i < parts.size(); i += 4)
				srcPoints.push_back(make_pair(atoi(parts[i].c_str()),
				Vector3D(atof(parts[i + 1].c_str()), atof(parts[i + 2].c_str()), atof(parts[i + 3].c_str()))));
		}
		else if (parts[0] == "dstVert")
		{
			for (int i = 1; i < parts.size(); ++i)
				dstVerts.push_back(atoi(parts[i].c_str()));
		}
		else if (parts[0] == "dstPoint")
		{
			for (int i = 1; i < parts.size(); i += 4)
				dstPoints.push_back(make_pair(atoi(parts[i].c_str()),
				Vector3D(atof(parts[i + 1].c_str()), atof(parts[i + 2].c_str()), atof(parts[i + 3].c_str()))));
		}
	}
	return true;
}

int main(int argc, char **argv)
{
	if (argc < 5)
	{
		cout << "USAGE: [.exe] [in.obj] [inputFile] [out.obj] [out.dist]" << endl;;
		return -1;
	}
	CMesh *mesh = new CMesh();
	if (!mesh->Load(argv[1]))
	{
		cout << "Cannot load mesh " << argv[1] << endl;
		return -2;
	}

	if (!LoadInputFile(argv[2]))
		return -3;

	ICH *ich = new ICH();

	ich->AssignMesh(mesh);
	for (int i = 0; i < srcVerts.size(); ++i)
		ich->AddSource(srcVerts[i]);
	for (int i = 0; i < srcPoints.size(); ++i)
		ich->AddSource(srcPoints[i].first, srcPoints[i].second);
	for (int i = 0; i < dstPoints.size(); ++i)
		ich->AddFacesKeptWindow(dstPoints[i].first);
	cout << "Executing ..." << endl;
	
	clock_t start = clock();
	ich->Execute(/*10000*/);
	clock_t end = clock();
	cout << "Time elapsed: " << (double)(end - start) / (double)CLOCKS_PER_SEC << endl;
	ich->OutputStatisticInfo();
// 	unsigned dstVert = 500, srcId;
// 	auto gp = ich->BuildGeodesicPathTo(dstVert, srcId);
// 	double totalLen = 0.0;
// 	Vector3D lastPoint = mesh->m_pVertex[dstVert].m_vPosition;
// 	for (auto iter = gp.begin(); iter != gp.end(); ++iter)
// 	{
// 		if (iter->isVertex)
// 		{
// 			cout << "Saddle vertex: " << iter->id << endl;
// 			totalLen += (lastPoint - mesh->m_pVertex[iter->id].m_vPosition).length();
// 			lastPoint = mesh->m_pVertex[iter->id].m_vPosition;
// 		}
// 		else
// 		{
// 			Vector3D p0 = mesh->m_pVertex[mesh->m_pEdge[iter->id].m_iVertex[0]].m_vPosition;
// 			Vector3D p1 = mesh->m_pVertex[mesh->m_pEdge[iter->id].m_iVertex[1]].m_vPosition;
// 			double l = mesh->m_pEdge[iter->id].m_length;
// 			Vector3D curPoint = (1.0 - iter->pos / l) * p0 + iter->pos / l * p1;
// 			cout << "Edge point: " << curPoint << endl;
// 			totalLen += (lastPoint - curPoint).length();
// 			lastPoint = curPoint;
// 		}
// 	}
// 	totalLen += (lastPoint - mesh->m_pVertex[srcId].m_vPosition).length();
// 	cout << "Accumulated length: " << totalLen << endl;
// 	cout << "Stored length: " << ich->GetDistanceTo(dstVert) << endl;

	// outputing
	cout << "Outputing ..." << endl;
	ofstream output(argv[argc - 2]);
	if (!output)
	{
		cout << "Cannot open output file " << argv[argc - 2] << endl;
		return -3;
	}
	output << "mtllib texture.mtl" << endl;
	for (int i = 0; i < mesh->m_nVertex; ++i)
		output << "v " << mesh->m_pVertex[i].m_vPosition << endl;
	unsigned errCnt = 0;
	double maxDist = 0.0;
	for (int i = 0; i < mesh->m_nVertex; ++i)
	{
		double dist = ich->GetDistanceTo(i);
		if (dist == DBL_MAX) continue;
		maxDist = max(maxDist, dist);
	}
	for (int i = 0; i < mesh->m_nVertex; ++i)
	{
		double dist = ich->GetDistanceTo(i);
		if (dist == DBL_MAX)
		{
			++errCnt;
			dist = 0.0;
		}
		output << "vt " << dist/maxDist << " " << dist/maxDist << endl;
	}
	for (int i = 0; i < mesh->m_nFace; ++i)
		output << "f " << mesh->m_pFace[i].m_piVertex[0] + 1 << "/" << mesh->m_pFace[i].m_piVertex[0] + 1 << " "
		<< mesh->m_pFace[i].m_piVertex[1] + 1 << "/" << mesh->m_pFace[i].m_piVertex[1] + 1 << " "
		<< mesh->m_pFace[i].m_piVertex[2] + 1 << "/" << mesh->m_pFace[i].m_piVertex[2] + 1 << endl;
	output.close();

	output.open(argv[argc - 1]);
	for (int i = 0; i < mesh->m_nVertex; ++i)
	{
		double dist = ich->GetDistanceTo(i);
		output << dist << endl;
	}
	output.close();

	cout << "Error cnt: " << errCnt << endl;

	for (int i = 0; i < dstVerts.size(); ++i)
	{
		int saddleCnt = 0;

		char outputFileName[255];
		sprintf_s(outputFileName, "%s.pathTo%d.vor", argv[1], dstVerts[i]);
		output.open(outputFileName);

		output << "0" << endl << "0" << endl;
		unsigned curSrcId = -1;
		auto gp = ich->BuildGeodesicPathTo(dstVerts[i], curSrcId);
		Vector3D curPoint = mesh->m_pVertex[dstVerts[i]].m_vPosition;
		Vector3D nextPoint;
		for (auto iter = gp.begin(); iter != gp.end(); ++iter)
		{
			if (iter->isVertex)
			{
				nextPoint = mesh->m_pVertex[iter->id].m_vPosition;
				++saddleCnt;
			}
			else
			{
				Vector3D p0 = mesh->m_pVertex[mesh->m_pEdge[iter->id].m_iVertex[0]].m_vPosition;
				Vector3D p1 = mesh->m_pVertex[mesh->m_pEdge[iter->id].m_iVertex[1]].m_vPosition;
				double l = mesh->m_pEdge[iter->id].m_length;
				nextPoint = (1.0 - iter->pos / l) * p0 + iter->pos / l * p1;
			}
			output << "face: 0 " << curPoint << " " << nextPoint << endl;
			curPoint = nextPoint;
		}
		nextPoint = curSrcId < mesh->m_nVertex ? mesh->m_pVertex[curSrcId].m_vPosition : srcPoints[curSrcId - mesh->m_nVertex].second;
		output << "face: 0 " << curPoint << " " << nextPoint << endl;
		output.close();
		cout << "Path to vert " << i << " passes " << saddleCnt << " saddle vertices." << endl;
	}

	for (int i = 0; i < dstPoints.size(); ++i)
	{
		int saddleCnt = 0;

		char outputFileName[255];
		sprintf_s(outputFileName, "%s.pathTo%d.vor", argv[1], mesh->m_nVertex + i);
		output.open(outputFileName);

		output << "0" << endl << "0" << endl;
		unsigned curSrcId = -1;
		auto gp = ich->BuildGeodesicPathTo(dstPoints[i].first, dstPoints[i].second, curSrcId);
		Vector3D curPoint = dstPoints[i].second;
		Vector3D nextPoint;
		for (auto iter = gp.begin(); iter != gp.end(); ++iter)
		{
			if (iter->isVertex)
			{
				nextPoint = mesh->m_pVertex[iter->id].m_vPosition;
				++saddleCnt;
			}
			else
			{
				Vector3D p0 = mesh->m_pVertex[mesh->m_pEdge[iter->id].m_iVertex[0]].m_vPosition;
				Vector3D p1 = mesh->m_pVertex[mesh->m_pEdge[iter->id].m_iVertex[1]].m_vPosition;
				double l = mesh->m_pEdge[iter->id].m_length;
				nextPoint = (1.0 - iter->pos / l) * p0 + iter->pos / l * p1;
			}
			output << "face: 0 " << curPoint << " " << nextPoint << endl;
			curPoint = nextPoint;
		}
		nextPoint = curSrcId < mesh->m_nVertex ? mesh->m_pVertex[curSrcId].m_vPosition : srcPoints[curSrcId - mesh->m_nVertex].second;
		output << "face: 0 " << curPoint << " " << nextPoint << endl;
		output.close();
		cout << "Path to point " << i << " passes " << saddleCnt << " saddle vertices." << endl;
	}
	return 0;
}