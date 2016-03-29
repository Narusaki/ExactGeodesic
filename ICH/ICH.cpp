// ICH.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "ICH.h"
#include <iostream>
#include <list>

using namespace std;

// This is the constructor of a class that has been exported.
// see ICH.h for the class definition
ICH::ICH()
{
	numOfWinGen = 0;
	maxWinQSize = 0;
	maxPseudoQSize = 0;
	totalCalcVertNum = 0;
	geodRadius = DBL_MAX;
	geodRadiusReached = false;
	return;
}

ICH::~ICH()
{
	return;
}

void ICH::Clear()
{
	while (!winQ.empty()) winQ.pop();
	while (!pseudoSrcQ.empty()) pseudoSrcQ.pop();

	for (int i = 0; i < mesh->m_nEdge; ++i)
	{
		splitInfos[i].dist = DBL_MAX;
		splitInfos[i].x = DBL_MAX;
		splitInfos[i].pseudoSrcId = -1;
	}
	for (int i = 0; i < mesh->m_nVertex; ++i)
	{
		vertInfos[i].birthTime = -1;
		vertInfos[i].dist = DBL_MAX;
		vertInfos[i].enterEdge = -1;
		vertInfos[i].isSource = false;
		vertInfos[i].pseudoSrcId = -1;
		vertInfos[i].srcId = -1;
	}
	sourceVerts.clear();
	sourcePoints.clear();

	storedWindows.clear();
	keptFaces.clear();

	numOfWinGen = 0;
	maxWinQSize = 0;
	maxPseudoQSize = 0;
	totalCalcVertNum = 0;
	geodRadius = DBL_MAX;
	geodRadiusReached = false;
}

void ICH::AssignMesh(CMesh *mesh_)
{
	mesh = mesh_;
	splitInfos.clear(); vertInfos.clear();
	splitInfos.resize(mesh->m_nEdge); vertInfos.resize(mesh->m_nVertex);
}

void ICH::AddSource(unsigned vertId)
{
	sourceVerts.push_back(vertId);
}

void ICH::AddSource(unsigned faceId, Vector3D pos)
{
	sourcePoints.push_back(make_pair(faceId, pos));
}

void ICH::AddFacesKeptWindow(unsigned faceId)
{
	keptFaces.push_back(faceId);
}

void ICH::SetMaxGeodRadius(double geodRadius_)
{
	geodRadius = geodRadius_;
}

void ICH::Execute(int totalCalcVertNum_)
{
	// Initialize
	Initialize();

	while (!winQ.empty() || !pseudoSrcQ.empty())
	{
		// Get valid window (for window whose pseudoSrcBirthTime is not equal (which means smaller/older) than
		// the current one, it must be an old window, which can be safely skipped)
		/*cout << "\r" << winQ.size() << " " << pseudoSrcQ.size();*/
		maxWinQSize = max(maxWinQSize, winQ.size());
		maxPseudoQSize = max(maxPseudoQSize, pseudoSrcQ.size());

		while (!winQ.empty() && winQ.top().pseudoSrcId < mesh->m_nVertex && 
			winQ.top().pseudoSrcBirthTime != vertInfos[winQ.top().pseudoSrcId].birthTime)
			winQ.pop();

		while (!pseudoSrcQ.empty() && winQ.top().pseudoSrcId < mesh->m_nVertex && 
			pseudoSrcQ.top().pseudoBirthTime != vertInfos[pseudoSrcQ.top().vertID].birthTime)
			pseudoSrcQ.pop();

		if (!winQ.empty() && (pseudoSrcQ.empty() || winQ.top().minDist < pseudoSrcQ.top().dist))
		{
			Window win = winQ.top(); winQ.pop();
			if (win.level > mesh->m_nFace) continue;
			// save windows for arbitrary dst geodesic construction
			unsigned twinEdge = mesh->m_pEdge[win.edgeID].m_iTwinEdge;
			if (twinEdge != -1)
				if (find(keptFaces.begin(), keptFaces.end(), mesh->m_pEdge[twinEdge].m_iFace) != keptFaces.end())
					storedWindows.push_back(win);
			PropagateWindow(win);
		}
		else if (!pseudoSrcQ.empty() && (winQ.empty() || winQ.top().minDist >= pseudoSrcQ.top().dist))
		{
			PseudoWindow pseudoWin = pseudoSrcQ.top(); pseudoSrcQ.pop();
			if (pseudoWin.level >= mesh->m_nFace) continue;
			GenSubWinsForPseudoSrc(pseudoWin);
		}

		if (totalCalcVertNum_ != -1 && totalCalcVertNum >= totalCalcVertNum_)
			break;
		if (geodRadiusReached)
			break;
	}
}

void ICH::OutputStatisticInfo()
{
	cout << "Total generated window number: " << numOfWinGen << endl;
	cout << "Max windows queue size: " << maxWinQSize << endl;
	cout << "Max pseudo-source queue size: " << maxPseudoQSize << endl;
	cout << "# of windows kept for arbitrary dst: " << storedWindows.size() << endl;
}

list<ICH::GeodesicKeyPoint> ICH::BuildGeodesicPathTo(unsigned faceId, Vector3D pos, unsigned &srcId)
{
	// find the window provide the nearest distance
	double minDist = DBL_MAX;
	Window minWin; double xInter; Vector2D pos2D;
	unsigned dstVert = -1;
	bool throughAWindow = true;

	// traverse the surrounded windows
	for (auto iter = storedWindows.begin(); iter != storedWindows.end(); ++iter)
	{
		unsigned twinEdge = mesh->m_pEdge[iter->edgeID].m_iTwinEdge;
		if (twinEdge == -1) continue;
		if (mesh->m_pEdge[twinEdge].m_iFace != faceId) continue;

		unsigned e0 = twinEdge;
		unsigned e1 = mesh->m_pEdge[e0].m_iNextEdge;
		unsigned e2 = mesh->m_pEdge[e1].m_iNextEdge;

		double l0 = mesh->m_pEdge[e0].m_length;
		double l1 = mesh->m_pEdge[e1].m_length;
		double l2 = mesh->m_pEdge[e2].m_length;

		unsigned v0 = mesh->m_pEdge[e0].m_iVertex[1];
		unsigned v1 = mesh->m_pEdge[e0].m_iVertex[0];
		unsigned v2 = mesh->m_pEdge[e1].m_iVertex[1];

		Vector2D p0(0.0, 0.0), p1(l0, 0.0), p2;
		p2.x = (l1*l1 + l0*l0 - l2*l2) / (2.0*l0);
		p2.y = -sqrt(fabs(l1*l1 - p2.x*p2.x));

		// window's pseudo source's 2D planar coordinate
		Vector2D src2D = iter->FlatenedSrc();

		// dst point's centroid coordinates
		double a = (pos - mesh->m_pVertex[v0].m_vPosition).length();
		double b = (pos - mesh->m_pVertex[v1].m_vPosition).length();
		double c = (pos - mesh->m_pVertex[v2].m_vPosition).length();

		double s0 = (b + c + l2) / 2.0;
		double s1 = (a + c + l1) / 2.0;
		double s2 = (a + b + l0) / 2.0;

		s0 = sqrt(fabs(s0 * (s0 - b) * (s0 - c) * (s0 - l2)));
		s1 = sqrt(fabs(s1 * (s1 - a) * (s1 - c) * (s1 - l1)));
		s2 = sqrt(fabs(s2 * (s2 - a) * (s2 - b) * (s2 - l0)));

		double w0 = s0 / (s0 + s1 + s2);
		double w1 = s1 / (s0 + s1 + s2);
		double w2 = s2 / (s0 + s1 + s2);

		Vector2D curPos2D = w0 * p0 + w1 * p1 + w2 * p2;

		// calculate the shortest distance
		double curXInter = src2D.x - (curPos2D.x - src2D.x) / (curPos2D.y - src2D.y) * src2D.y;
		if (curXInter <= 0.0 || curXInter >= l0) continue;
		double curMinDist = DBL_MAX;
		if (curXInter > iter->b0 && curXInter < iter->b1)
			curMinDist = (curPos2D - src2D).length() + iter->pseudoSrcDist;
		else if (curXInter <= iter->b0)
			curMinDist = (curPos2D - Vector2D(iter->b0, 0.0)).length() + iter->d0 + iter->pseudoSrcDist;
		else
			curMinDist = (curPos2D - Vector2D(iter->b1, 0.0)).length() + iter->d1 + iter->pseudoSrcDist;

		if (curMinDist < minDist)
		{
			minDist = curMinDist;
			minWin = *iter;
			xInter = curXInter;
			pos2D = curPos2D;
		}
	}

	// traverse the surrounded vertices
	for (int i = 0; i < 3; ++i)
	{
		unsigned opVert = mesh->m_pFace[faceId].m_piVertex[i];
		if (mesh->m_pAngles[opVert] < 2.0 * PI) continue;

		double curDist = (pos - mesh->m_pVertex[opVert].m_vPosition).length() + vertInfos[opVert].dist;
		if (curDist < minDist)
		{
			throughAWindow = false;
			dstVert = opVert;
			minDist = curDist;
		}
	}

	if (!throughAWindow)
	{
		auto path = BuildGeodesicPathTo(dstVert, srcId);
		GeodesicKeyPoint gkp;
		gkp.isVertex = true; gkp.id = dstVert;
		path.push_front(gkp);
		return path;
	}
	else
	{
		// next key point is on an edge
		list< GeodesicKeyPoint > path;
		GeodesicKeyPoint gkp;
		gkp.isVertex = false;
		gkp.id = mesh->m_pEdge[minWin.edgeID].m_iTwinEdge; 
		gkp.pos = mesh->m_pEdge[gkp.id].m_length - xInter;
		path.push_back(gkp);

		unsigned enterEdge = gkp.id;
		unsigned opVert = mesh->m_pEdge[gkp.id].m_iTwinEdge;
		opVert = mesh->m_pEdge[mesh->m_pEdge[opVert].m_iNextEdge].m_iVertex[1];
		double l0 = mesh->m_pEdge[gkp.id].m_length;
		double l1 = mesh->m_pEdge[mesh->m_pEdge[gkp.id].m_iNextEdge].m_length;
		double l2 = mesh->m_pEdge[mesh->m_pEdge[mesh->m_pEdge[gkp.id].m_iNextEdge].m_iNextEdge].m_length;

		Vector2D lastPoint = pos2D, curPoint;
		curPoint.x = l0 - gkp.pos; curPoint.y = 0.0;

		while (minWin.pseudoSrcId < mesh->m_nVertex && opVert != minWin.pseudoSrcId ||
			minWin.pseudoSrcId >= mesh->m_nVertex &&
			mesh->m_pEdge[mesh->m_pEdge[gkp.id].m_iTwinEdge].m_iFace != sourcePoints[minWin.pseudoSrcId - mesh->m_nVertex].first)
		{
			// trace back
			unsigned e0 = mesh->m_pEdge[gkp.id].m_iTwinEdge;
			unsigned e1 = mesh->m_pEdge[e0].m_iNextEdge;
			unsigned e2 = mesh->m_pEdge[e1].m_iNextEdge;
			double l0 = mesh->m_pEdge[e0].m_length;
			double l1 = mesh->m_pEdge[e1].m_length;
			double l2 = mesh->m_pEdge[e2].m_length;

			Vector2D opVert2D;
			opVert2D.x = (l0*l0 + l2*l2 - l1*l1) / (2.0*l0);
			opVert2D.y = sqrt(fabs(l2*l2 - opVert2D.x*opVert2D.x));

			if (toLeft(opVert2D, lastPoint, curPoint))
			{
				Vector2D p0, p1;
				p0.x = (l2*l2 + l1*l1 - l0*l0) / (2.0*l1);
				p0.y = -sqrt(fabs(l2*l2 - p0.x*p0.x));
				p1.x = l1; p1.y = 0.0;
				Vector2D newlastPoint = gkp.pos / l0 * p0 + (1.0 - gkp.pos / l0) * p1;

				gkp.pos = Intersect(lastPoint, curPoint, Vector2D(l0, 0.0), opVert2D);
				gkp.pos = (1.0 - gkp.pos) * l1;
				gkp.id = e1;
				curPoint.x = l1 - gkp.pos; curPoint.y = 0.0;
				lastPoint = newlastPoint;
			}
			else
			{
				Vector2D p0, p1;
				p0.x = 0.0; p0.y = 0.0;
				p1.x = (l2*l2 + l0*l0 - l1*l1) / (2.0*l2);
				p1.y = -sqrt(fabs(l0*l0 - p1.x*p1.x));
				Vector2D newlastPoint = gkp.pos / l0 * p0 + (1.0 - gkp.pos / l0) * p1;

				gkp.pos = Intersect(lastPoint, curPoint, opVert2D, Vector2D(0.0, 0.0));
				gkp.pos = (1.0 - gkp.pos) * l2;
				gkp.id = e2;
				curPoint.x = l2 - gkp.pos; curPoint.y = 0.0;
				lastPoint = newlastPoint;
			}
			path.push_back(gkp);

			opVert = mesh->m_pEdge[gkp.id].m_iTwinEdge;
			opVert = mesh->m_pEdge[mesh->m_pEdge[opVert].m_iNextEdge].m_iVertex[1];
		}

		if (minWin.pseudoSrcId >= mesh->m_nVertex) {
			dstVert = minWin.pseudoSrcId;
			srcId = dstVert;
		}
		else if (!vertInfos[opVert].isSource)
		{
			gkp.isVertex = true;
			gkp.id = opVert;
			path.push_back(gkp);
			dstVert = opVert;

			auto subPath = BuildGeodesicPathTo(opVert, srcId);
			path.insert(path.end(), subPath.begin(), subPath.end());
		}
		srcId = minWin.srcID;
		return path;
	}
}

list<ICH::GeodesicKeyPoint> ICH::BuildGeodesicPathTo(unsigned vertId, unsigned &srcId)
{
	// TODO: build geodesic path from vertex vertId to source
	list < GeodesicKeyPoint > path;
	unsigned curVert = vertId;
	GeodesicKeyPoint gkp;
	while (!vertInfos[curVert].isSource)
	{
		unsigned enterEdge = vertInfos[curVert].enterEdge;
		if (enterEdge == -1)
		{
			// trace back to an arbitrary point
			double planarDist = DBL_MAX;
			for (int i = 0; i < sourcePoints.size(); ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					if (mesh->m_pFace[sourcePoints[i].first].m_piVertex[j] != curVert) continue;
					double curPlanarDist = (sourcePoints[i].second - mesh->m_pVertex[curVert].m_vPosition).length();
					planarDist = min(planarDist, curPlanarDist);
					srcId = mesh->m_nVertex + i;
					break;
				}
			}
			return path;
		}
		else if (mesh->m_pEdge[enterEdge].m_iVertex[0] == curVert)
		{
			// next key point is still a vertex
			unsigned nextVert = mesh->m_pEdge[enterEdge].m_iVertex[1];
			if (!vertInfos[nextVert].isSource)
			{
				gkp.isVertex = true;
				gkp.id = nextVert;
				path.push_back(gkp);
			}
			curVert = nextVert;
		}
		else
		{
			// next key point is on an edge
			gkp.isVertex = false;
			gkp.id = enterEdge; gkp.pos = splitInfos[enterEdge].x;
			path.push_back(gkp);

			unsigned opVert = mesh->m_pEdge[gkp.id].m_iTwinEdge;
			opVert = mesh->m_pEdge[mesh->m_pEdge[opVert].m_iNextEdge].m_iVertex[1];
			double l0 = mesh->m_pEdge[gkp.id].m_length;
			double l1 = mesh->m_pEdge[mesh->m_pEdge[gkp.id].m_iNextEdge].m_length;
			double l2 = mesh->m_pEdge[mesh->m_pEdge[mesh->m_pEdge[gkp.id].m_iNextEdge].m_iNextEdge].m_length;
			
			Vector2D lastPoint, curPoint;
			lastPoint.x = (l1*l1 + l0*l0 - l2*l2) / (2.0*l0);
			lastPoint.y = -sqrt(fabs(l1*l1 - lastPoint.x*lastPoint.x));
			curPoint.x = l0 - gkp.pos; curPoint.y = 0.0;

			while (splitInfos[enterEdge].pseudoSrcId < mesh->m_nVertex && 
				opVert != splitInfos[enterEdge].pseudoSrcId || 
				splitInfos[enterEdge].pseudoSrcId >= mesh->m_nVertex && 
				mesh->m_pEdge[mesh->m_pEdge[gkp.id].m_iTwinEdge].m_iFace != sourcePoints[splitInfos[enterEdge].pseudoSrcId - mesh->m_nVertex].first)
			{
				// trace back
				unsigned e0 = mesh->m_pEdge[gkp.id].m_iTwinEdge;
				unsigned e1 = mesh->m_pEdge[e0].m_iNextEdge;
				unsigned e2 = mesh->m_pEdge[e1].m_iNextEdge;
				double l0 = mesh->m_pEdge[e0].m_length;
				double l1 = mesh->m_pEdge[e1].m_length;
				double l2 = mesh->m_pEdge[e2].m_length;

				Vector2D opVert2D;
				opVert2D.x = (l0*l0 + l2*l2 - l1*l1) / (2.0*l0);
				opVert2D.y = sqrt(fabs(l2*l2 - opVert2D.x*opVert2D.x));

				if (toLeft(opVert2D, lastPoint, curPoint))
				{
					Vector2D p0, p1;
					p0.x = (l2*l2 + l1*l1 - l0*l0) / (2.0*l1);
					p0.y = -sqrt(fabs(l2*l2 - p0.x*p0.x));
					p1.x = l1; p1.y = 0.0;
					Vector2D newlastPoint = gkp.pos / l0 * p0 + (1.0 - gkp.pos / l0) * p1;

					gkp.pos = Intersect(lastPoint, curPoint, Vector2D(l0, 0.0), opVert2D);
					gkp.pos = (1.0 - gkp.pos) * l1;
					gkp.id = e1;
					curPoint.x = l1 - gkp.pos; curPoint.y = 0.0;
					lastPoint = newlastPoint;
				}
				else
				{
					Vector2D p0, p1;
					p0.x = 0.0; p0.y = 0.0;
					p1.x = (l2*l2 + l0*l0 - l1*l1) / (2.0*l2);
					p1.y = -sqrt(fabs(l0*l0 - p1.x*p1.x));
					Vector2D newlastPoint = gkp.pos / l0 * p0 + (1.0 - gkp.pos / l0) * p1;

					gkp.pos = Intersect(lastPoint, curPoint, opVert2D, Vector2D(0.0, 0.0));
					gkp.pos = (1.0 - gkp.pos) * l2;
					gkp.id = e2;
					curPoint.x = l2 - gkp.pos; curPoint.y = 0.0;
					lastPoint = newlastPoint;
				}
				path.push_back(gkp);

				opVert = mesh->m_pEdge[gkp.id].m_iTwinEdge;
				opVert = mesh->m_pEdge[mesh->m_pEdge[opVert].m_iNextEdge].m_iVertex[1];
			}

			if (splitInfos[enterEdge].pseudoSrcId >= mesh->m_nVertex) {
				curVert = splitInfos[enterEdge].pseudoSrcId;
				break;
			}
			if (!vertInfos[opVert].isSource)
			{
				gkp.isVertex = true;
				gkp.id = opVert;
				path.push_back(gkp);
			}
			curVert = opVert;
		}
	}
	srcId = curVert;
	return path;
}

double ICH::GetDistanceTo(unsigned vertId)
{
	return vertInfos[vertId].dist;
}

void ICH::Initialize()
{
	for (int i = 0; i < sourceVerts.size(); ++i)
	{
		unsigned srcId = sourceVerts[i];
		for (int j = 0; j < mesh->m_pVertex[srcId].m_nValence; ++j)
		{
			unsigned opEdge = mesh->m_pEdge[mesh->m_pVertex[srcId].m_piEdge[j]].m_iNextEdge;
			Window win;
			win.edgeID = opEdge;
			win.b0 = 0.0; win.b1 = mesh->m_pEdge[opEdge].m_length;
			win.d0 = mesh->m_pEdge[mesh->m_pVertex[srcId].m_piEdge[j]].m_length;
			win.d1 = mesh->m_pEdge[mesh->m_pEdge[opEdge].m_iNextEdge].m_length;
			win.pseudoSrcDist = 0.0; win.calcMinDist();
			win.srcID = srcId; win.pseudoSrcId = srcId;
			win.pseudoSrcBirthTime = 0;		
			win.level = 0;
			winQ.push(win);
			++numOfWinGen;

			unsigned opVert = mesh->m_pEdge[mesh->m_pVertex[srcId].m_piEdge[j]].m_iVertex[1];
			vertInfos[opVert].birthTime = 0;
			vertInfos[opVert].dist = mesh->m_pEdge[mesh->m_pVertex[srcId].m_piEdge[j]].m_length;
			vertInfos[opVert].enterEdge = mesh->m_pEdge[mesh->m_pVertex[srcId].m_piEdge[j]].m_iTwinEdge;
			vertInfos[opVert].srcId = srcId; vertInfos[opVert].pseudoSrcId = srcId;

			if (mesh->m_pAngles[opVert] < 2.0 * PI) continue;

			PseudoWindow pseudoWin;
			pseudoWin.vertID = opVert; pseudoWin.dist = mesh->m_pEdge[mesh->m_pVertex[srcId].m_piEdge[j]].m_length;
			pseudoWin.srcId = srcId; pseudoWin.pseudoSrcId = srcId;
			pseudoWin.pseudoBirthTime = vertInfos[opVert].birthTime;
			pseudoWin.level = 0;
			pseudoSrcQ.push(pseudoWin);
		}
		vertInfos[srcId].birthTime = 0;
		vertInfos[srcId].dist = 0.0;
		vertInfos[srcId].enterEdge = -1;
		vertInfos[srcId].isSource = true;
		vertInfos[srcId].srcId = srcId; vertInfos[srcId].pseudoSrcId = srcId;
	}

	for (int i = 0; i < sourcePoints.size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			unsigned opEdge = mesh->m_pFace[sourcePoints[i].first].m_piEdge[j];
			Window win;
			win.edgeID = opEdge;
			win.b0 = 0.0; win.b1 = mesh->m_pEdge[opEdge].m_length;
			win.d0 = (sourcePoints[i].second - mesh->m_pVertex[mesh->m_pEdge[opEdge].m_iVertex[0]].m_vPosition).length();
			win.d1 = (sourcePoints[i].second - mesh->m_pVertex[mesh->m_pEdge[opEdge].m_iVertex[1]].m_vPosition).length();
			win.pseudoSrcDist = 0.0; win.calcMinDist();
			win.srcID = mesh->m_nVertex + i; win.pseudoSrcId = win.srcID;
			win.pseudoSrcBirthTime = 0; win.level = 0;
			winQ.push(win);

			unsigned opVert = mesh->m_pEdge[opEdge].m_iVertex[0];
			vertInfos[opVert].birthTime = 0;
			vertInfos[opVert].dist = (sourcePoints[i].second - mesh->m_pVertex[opVert].m_vPosition).length();
			vertInfos[opVert].enterEdge = -1;
			vertInfos[opVert].srcId = mesh->m_nVertex + i; vertInfos[opVert].pseudoSrcId = mesh->m_nVertex + i;

			if (mesh->m_pAngles[opVert] < 2.0 * PI) continue;

			PseudoWindow pseudoWin;
			pseudoWin.vertID = opVert; 
			pseudoWin.dist = (mesh->m_pVertex[opVert].m_vPosition - sourcePoints[i].second).length();
			pseudoWin.srcId = win.srcID; pseudoWin.pseudoSrcId = win.srcID;
			pseudoWin.pseudoBirthTime = vertInfos[opVert].birthTime;
			pseudoWin.level = 0;
			pseudoSrcQ.push(pseudoWin);
		}
	}
}

void ICH::PropagateWindow(const Window &win)
{
	unsigned e0 = mesh->m_pEdge[win.edgeID].m_iTwinEdge;
	if (e0 == -1) return;
	unsigned e1 = mesh->m_pEdge[e0].m_iNextEdge;
	unsigned e2 = mesh->m_pEdge[e1].m_iNextEdge;
	unsigned opVert = mesh->m_pEdge[e1].m_iVertex[1];

	Vector2D src2D = win.FlatenedSrc();

	Vector2D left(win.b0, 0.0), right(win.b1, 0.0);
	double l0 = mesh->m_pEdge[e0].m_length;
	double l1 = mesh->m_pEdge[e1].m_length;
	double l2 = mesh->m_pEdge[e2].m_length;
	Vector2D v0(0.0, 0.0), v1(l0, 0.0), v2;
	v2.x = (l1*l1 + l0*l0 - l2*l2) / (2.0 * l0);
	v2.y = -sqrt(fabs(l1*l1 - v2.x*v2.x));

	double interX = v2.x - v2.y * (v2.x - src2D.x) / (v2.y - src2D.y);
	Window leftChildWin, rightChildWin;
	bool hasLeftChild = true, hasRightChild = true;
	// only generate right window
	if (interX <= left.x)
	{
		hasLeftChild = false;
		double t0 = Intersect(src2D, left, v2, v1);
		double t1 = Intersect(src2D, right, v2, v1);
		BuildWindow(win, e2, t0, t1, v2, v1, rightChildWin);
		if (!IsValidWindow(rightChildWin, false)) hasRightChild = false;
	}
	// only generate left window
	else if (interX >= right.x)
	{
		hasRightChild = false;
		double t0 = Intersect(src2D, left, v0, v2);
		double t1 = Intersect(src2D, right, v0, v2);
		BuildWindow(win, e1, t0, t1, v0, v2, leftChildWin);
		if (!IsValidWindow(leftChildWin, true)) hasLeftChild = false;
	}
	// generate both left and right window
	else
	{
		double directDist = (v2 - src2D).length();
		// ONE ANGLE, ONE SPLIT
		if (directDist + win.pseudoSrcDist > splitInfos[e0].dist && 
			(directDist + win.pseudoSrcDist) / splitInfos[e0].dist - 1.0 > RELATIVE_ERROR)
		{
			hasLeftChild = splitInfos[e0].x < interX;
			hasRightChild = !hasLeftChild;
			/*cout << "Filter 1 works..." << endl;*/
		}
		else
		{
			if (directDist + win.pseudoSrcDist < splitInfos[e0].dist)
			{
				splitInfos[e0].dist = directDist + win.pseudoSrcDist;
				splitInfos[e0].pseudoSrcId = win.pseudoSrcId;
				splitInfos[e0].srcId = win.srcID;
				splitInfos[e0].level = win.level;
				splitInfos[e0].x = l0 - interX;
			}

			if (directDist + win.pseudoSrcDist < vertInfos[opVert].dist)
			{
				if (vertInfos[opVert].dist == DBL_MAX)
					++totalCalcVertNum;
				if (directDist + win.pseudoSrcDist >= geodRadius)
					geodRadiusReached = true;

				++vertInfos[opVert].birthTime;
				vertInfos[opVert].dist = directDist + win.pseudoSrcDist;
				vertInfos[opVert].enterEdge = e0;
				vertInfos[opVert].srcId = win.srcID; vertInfos[opVert].pseudoSrcId = win.pseudoSrcId;
				if (mesh->m_pAngles[opVert] > 2.0 * PI)
				{
					PseudoWindow pseudoWin;
					pseudoWin.vertID = opVert; pseudoWin.dist = vertInfos[opVert].dist;
					pseudoWin.srcId = win.srcID; pseudoWin.pseudoSrcId = win.pseudoSrcId;
					pseudoWin.pseudoBirthTime = vertInfos[opVert].birthTime;
					pseudoWin.level = win.level + 1;
					pseudoSrcQ.push(pseudoWin);
				}
			}
		}
		if (hasLeftChild)
		{
			// left child window
			double t0 = Intersect(src2D, left, v0, v2);
			BuildWindow(win, e1, t0, 0.0, v0, v2, leftChildWin);
			if (!IsValidWindow(leftChildWin, true)) hasLeftChild = false;
		}
		if (hasRightChild)
		{
			// right child window
			double t1 = Intersect(src2D, right, v2, v1);
			BuildWindow(win, e2, 1.0, t1, v2, v1, rightChildWin);
			if (!IsValidWindow(rightChildWin, false)) hasRightChild = false;
		}
	}

	if (hasLeftChild)
	{
		++numOfWinGen;
		winQ.push(leftChildWin);
	}
	if (hasRightChild)
	{
		++numOfWinGen;
		winQ.push(rightChildWin);
	}

}

void ICH::GenSubWinsForPseudoSrc(const PseudoWindow &pseudoWin)
{
	unsigned startEdge, endEdge;
	if (vertInfos[pseudoWin.vertID].enterEdge == -1 && vertInfos[pseudoWin.vertID].birthTime != -1)
	{
		startEdge = mesh->m_pVertex[pseudoWin.vertID].m_piEdge[0];
		endEdge = startEdge;
	}
	else if (mesh->m_pEdge[vertInfos[pseudoWin.vertID].enterEdge].m_iVertex[0] == pseudoWin.vertID)
		GenSubWinsForPseudoSrcFromPseudoSrc(pseudoWin, startEdge, endEdge);
	else if (mesh->m_pEdge[mesh->m_pEdge[vertInfos[pseudoWin.vertID].enterEdge].m_iNextEdge].m_iVertex[1] == pseudoWin.vertID)
		GenSubWinsForPseudoSrcFromWindow(pseudoWin, startEdge, endEdge);
	else assert(false);

	// generate windows
	do
	{
		Window win;
		win.edgeID = mesh->m_pEdge[startEdge].m_iNextEdge;
		win.b0 = 0.0; win.b1 = mesh->m_pEdge[win.edgeID].m_length;
		win.d0 = mesh->m_pEdge[startEdge].m_length;
		win.d1 = mesh->m_pEdge[mesh->m_pEdge[win.edgeID].m_iNextEdge].m_length;
		win.pseudoSrcDist = pseudoWin.dist; win.calcMinDist();
		win.srcID = pseudoWin.srcId; win.pseudoSrcId = pseudoWin.vertID;
		win.pseudoSrcBirthTime = pseudoWin.pseudoBirthTime;
		win.level = pseudoWin.level + 1;
		winQ.push(win);
		++numOfWinGen;

		startEdge = mesh->m_pEdge[mesh->m_pEdge[mesh->m_pEdge[startEdge].m_iNextEdge].m_iNextEdge].m_iTwinEdge;
	} while (startEdge != endEdge);

	// generate adjacent pseudo sources
	for (int i = 0; i < mesh->m_pVertex[pseudoWin.vertID].m_nValence; ++i)
	{
		unsigned adjEdge = mesh->m_pVertex[pseudoWin.vertID].m_piEdge[i];
		unsigned opVert = mesh->m_pEdge[adjEdge].m_iVertex[1];
		if (mesh->m_pAngles[opVert] < 2.0 * PI) continue;
		if (vertInfos[opVert].dist < pseudoWin.dist + mesh->m_pEdge[adjEdge].m_length) continue;

		if (vertInfos[opVert].dist == DBL_MAX)
			++totalCalcVertNum;
		if (pseudoWin.dist + mesh->m_pEdge[adjEdge].m_length >= geodRadius)
			geodRadiusReached = true;

		vertInfos[opVert].dist = pseudoWin.dist + mesh->m_pEdge[adjEdge].m_length;
		++vertInfos[opVert].birthTime;
		vertInfos[opVert].enterEdge = mesh->m_pEdge[adjEdge].m_iTwinEdge;
		vertInfos[opVert].srcId = pseudoWin.srcId; vertInfos[opVert].pseudoSrcId = pseudoWin.pseudoSrcId;

		PseudoWindow childPseudoWin;
		childPseudoWin.vertID = opVert; childPseudoWin.dist = vertInfos[opVert].dist;
		childPseudoWin.srcId = pseudoWin.srcId; childPseudoWin.pseudoSrcId = pseudoWin.vertID;
		childPseudoWin.pseudoBirthTime = vertInfos[opVert].birthTime;
		childPseudoWin.level = pseudoWin.level;
		pseudoSrcQ.push(childPseudoWin);
	}
}

void ICH::GenSubWinsForPseudoSrcFromWindow(const PseudoWindow &pseudoWin, unsigned &startEdge, unsigned &endEdge)
{
	unsigned e0 = vertInfos[pseudoWin.vertID].enterEdge;
	unsigned e1 = mesh->m_pEdge[e0].m_iNextEdge;
	unsigned e2 = mesh->m_pEdge[e1].m_iNextEdge;

	double l0 = mesh->m_pEdge[e0].m_length;
	double l1 = mesh->m_pEdge[e1].m_length;
	double l2 = mesh->m_pEdge[e2].m_length;

	unsigned pseudoSrc = pseudoWin.vertID;
	Vector2D enterPoint;
	enterPoint.x = l0 - splitInfos[e0].x;
	enterPoint.y = 0.0;

	Vector2D v0(0.0, 0.0), v1(l0, 0.0), v2;
	v2.x = (l1*l1 + l0*l0 - l2*l2) / (2.0*l0);
	v2.y = -sqrt(fabs(l1*l1 - v2.x*v2.x));

	// TODO: generate windows using opVert as pseudo sources
	double angle0 = (enterPoint - v2) * (v0 - v2) / (enterPoint - v2).length() / l1;
	double angle1 = (enterPoint - v2) * (v1 - v2) / (enterPoint - v2).length() / l2;
	if (angle0 > 1.0) angle0 = 1.0; else if (angle0 < -1.0) angle0 = -1.0;
	if (angle1 > 1.0) angle1 = 1.0; else if (angle1 < -1.0) angle1 = -1.0;
	angle0 = acos(angle0); angle1 = acos(angle1);

	startEdge = -1, endEdge = -1;
	// traverse from left
	unsigned curEdge = mesh->m_pEdge[e1].m_iTwinEdge;
	while (angle0 < PI && curEdge != -1)
	{
		unsigned opEdge = mesh->m_pEdge[curEdge].m_iNextEdge;
		unsigned nextEdge = mesh->m_pEdge[opEdge].m_iNextEdge;
		double L0 = mesh->m_pEdge[curEdge].m_length, L1 = mesh->m_pEdge[nextEdge].m_length;
		double L2 = mesh->m_pEdge[opEdge].m_length;
		double curAngle = (L0*L0 + L1*L1 - L2*L2) / (2.0 * L0 * L1);
		if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
		curAngle = acos(curAngle);
		angle0 += curAngle;
		curEdge = mesh->m_pEdge[nextEdge].m_iTwinEdge;
	}
	if (curEdge != -1)
		startEdge = mesh->m_pEdge[mesh->m_pEdge[curEdge].m_iTwinEdge].m_iNextEdge;

	// traverse from right
	curEdge = mesh->m_pEdge[e2].m_iTwinEdge;
	while (angle1 < PI && curEdge != -1)
	{
		unsigned nextEdge = mesh->m_pEdge[curEdge].m_iNextEdge;
		unsigned opEdge = mesh->m_pEdge[nextEdge].m_iNextEdge;
		double L0 = mesh->m_pEdge[curEdge].m_length, L1 = mesh->m_pEdge[nextEdge].m_length;
		double L2 = mesh->m_pEdge[opEdge].m_length;
		double curAngle = (L0*L0 + L1*L1 - L2*L2) / (2.0 * L0 * L1);
		if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
		curAngle = acos(curAngle);
		angle1 += curAngle;
		curEdge = mesh->m_pEdge[nextEdge].m_iTwinEdge;
	}
	if (curEdge != -1)
	{
		endEdge = mesh->m_pEdge[mesh->m_pEdge[mesh->m_pEdge[curEdge].m_iTwinEdge].m_iNextEdge].m_iNextEdge;
		endEdge = mesh->m_pEdge[endEdge].m_iTwinEdge;
	}
}

void ICH::GenSubWinsForPseudoSrcFromPseudoSrc(const PseudoWindow &pseudoWin, unsigned &startEdge, unsigned &endEdge)
{
	unsigned pseudoSrc = pseudoWin.vertID;

	// TODO: generate windows using opVert as pseudo sources
	double angle0 = 0.0, angle1 = 0.0;

	startEdge = -1, endEdge = -1;
	// traverse from left
	unsigned curEdge = vertInfos[pseudoWin.vertID].enterEdge;
	while (angle0 < PI && curEdge != -1)
	{
		unsigned opEdge = mesh->m_pEdge[curEdge].m_iNextEdge;
		unsigned nextEdge = mesh->m_pEdge[opEdge].m_iNextEdge;
		double L0 = mesh->m_pEdge[curEdge].m_length, L1 = mesh->m_pEdge[nextEdge].m_length;
		double L2 = mesh->m_pEdge[opEdge].m_length;
		double curAngle = (L0*L0 + L1*L1 - L2*L2) / (2.0 * L0 * L1);
		if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
		curAngle = acos(curAngle);
		angle0 += curAngle;
		curEdge = mesh->m_pEdge[nextEdge].m_iTwinEdge;
	}
	if (curEdge != -1)
		startEdge = mesh->m_pEdge[mesh->m_pEdge[curEdge].m_iTwinEdge].m_iNextEdge;

	// traverse from right
	curEdge = mesh->m_pEdge[vertInfos[pseudoWin.vertID].enterEdge].m_iTwinEdge;
	while (angle1 < PI && curEdge != -1)
	{
		unsigned nextEdge = mesh->m_pEdge[curEdge].m_iNextEdge;
		unsigned opEdge = mesh->m_pEdge[nextEdge].m_iNextEdge;
		double L0 = mesh->m_pEdge[curEdge].m_length, L1 = mesh->m_pEdge[nextEdge].m_length;
		double L2 = mesh->m_pEdge[opEdge].m_length;
		double curAngle = (L0*L0 + L1*L1 - L2*L2) / (2.0 * L0 * L1);
		if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
		curAngle = acos(curAngle);
		angle1 += curAngle;
		curEdge = mesh->m_pEdge[nextEdge].m_iTwinEdge;
	}
	if (curEdge != -1)
	{
		endEdge = mesh->m_pEdge[mesh->m_pEdge[mesh->m_pEdge[curEdge].m_iTwinEdge].m_iNextEdge].m_iNextEdge;
		endEdge = mesh->m_pEdge[endEdge].m_iTwinEdge;
	}
}

bool ICH::IsValidWindow(const Window &win, bool isLeftChild)
{
	if (win.b1 <= win.b0) return false;
	// apply ICH's filter
	unsigned v1 = mesh->m_pEdge[win.edgeID].m_iVertex[0];
	unsigned v2 = mesh->m_pEdge[win.edgeID].m_iVertex[1];
	unsigned v3 = mesh->m_pEdge[mesh->m_pEdge[win.edgeID].m_iNextEdge].m_iVertex[1];
	double l0 = mesh->m_pEdge[win.edgeID].m_length;
	double l1 = mesh->m_pEdge[mesh->m_pEdge[win.edgeID].m_iNextEdge].m_length;
	double l2 = mesh->m_pEdge[mesh->m_pEdge[mesh->m_pEdge[win.edgeID].m_iNextEdge].m_iNextEdge].m_length;
	Vector2D p1(0.0, 0.0), p2(l0, 0.0), p3;
	p3.x = (l2*l2 + l0*l0 - l1*l1) / (2.0 * l0);
	p3.y = sqrt(fabs(l2*l2 - p3.x*p3.x));

	Vector2D A(win.b0, 0.0), B(win.b1, 0.0);
	Vector2D src2D = win.FlatenedSrc();
	

	if (win.pseudoSrcDist + (src2D - B).length() > vertInfos[v1].dist + win.b1 && 
		(win.pseudoSrcDist + (src2D - B).length()) / (vertInfos[v1].dist + win.b1) - 1.0 > 0.0)
	{
		/*cout << "Filter 2 works..." << endl;*/
		return false;
	}
	if (win.pseudoSrcDist + (src2D - A).length() > vertInfos[v2].dist + l0 - win.b0 && 
		(win.pseudoSrcDist + (src2D - A).length()) / (vertInfos[v2].dist + l0 - win.b0) - 1.0 > 0.0)
	{
		/*cout << "Filter 2 works..." << endl;*/
		return false;
	}
	if (isLeftChild)
	{
		if (win.pseudoSrcDist + (src2D - A).length() > vertInfos[v3].dist + (p3 - A).length() && 
			(win.pseudoSrcDist + (src2D - A).length()) / (vertInfos[v3].dist + (p3 - A).length()) - 1.0 > 0.0)
		{
			/*cout << "Filter 2 works..." << endl;*/
			return false;
		}
	}
	else
	{
		if (win.pseudoSrcDist + (src2D - B).length() > vertInfos[v3].dist + (p3 - B).length() && 
			(win.pseudoSrcDist + (src2D - B).length()) / (vertInfos[v3].dist + (p3 - B).length()) - 1.0 > RELATIVE_ERROR)
		{
			/*cout << "Filter 2 works..." << endl;*/
			return false;
		}
	}
	return true;
}

void ICH::BuildWindow(const Window &fatherWin, 
	unsigned edge, 
	double t0, double t1, 
	const Vector2D &v0, const Vector2D &v1, 
	Window &win)
{
	Vector2D src2D = fatherWin.FlatenedSrc();
	win.edgeID = edge;
	win.b0 = (1 - t0) * mesh->m_pEdge[edge].m_length; win.b1 = (1 - t1) * mesh->m_pEdge[edge].m_length;
	win.d0 = (src2D - (t0 * v0 + (1 - t0)*v1)).length();
	win.d1 = (src2D - (t1 * v0 + (1 - t1)*v1)).length();
	win.pseudoSrcDist = fatherWin.pseudoSrcDist;
	win.calcMinDist();
	win.srcID = fatherWin.srcID; win.pseudoSrcId = fatherWin.pseudoSrcId;
	win.pseudoSrcBirthTime = fatherWin.pseudoSrcBirthTime;
	win.level = fatherWin.level + 1;
}

double ICH::Intersect(const Vector2D &v0, const Vector2D &v1, const Vector2D &p0, const Vector2D &p1)
{
	double a00 = p0.x - p1.x, a01 = v1.x - v0.x, b0 = v1.x - p1.x;
	double a10 = p0.y - p1.y, a11 = v1.y - v0.y, b1 = v1.y - p1.y;
	return (b0*a11 - b1*a01) / (a00*a11 - a10*a01);
}