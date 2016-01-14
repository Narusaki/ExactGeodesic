// ICH.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "ICH.h"
#include <iostream>

using namespace std;

// This is the constructor of a class that has been exported.
// see ICH.h for the class definition
ICH::ICH()
{
	return;
}

ICH::~ICH()
{
	return;
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

void ICH::Execute()
{
	// Initialize
	Initialize();

	while (!winQ.empty() || !pseudoSrcQ.empty())
	{
		// Get valid window (for window whose pseudoSrcBirthTime is not equal (which means smaller/older) than
		// the current one, it must be an old window, which can be safely skipped)
		/*cout << "\r" << winQ.size() << " " << pseudoSrcQ.size();*/
		while (!winQ.empty() &&
			winQ.top().pseudoSrcBirthTime != vertInfos[winQ.top().pseudoSrcId].birthTime)
			winQ.pop();

		while (!pseudoSrcQ.empty() &&
			pseudoSrcQ.top().pseudoBirthTime != vertInfos[pseudoSrcQ.top().vertID].birthTime)
			pseudoSrcQ.pop();

		if (!winQ.empty() && (pseudoSrcQ.empty() || winQ.top().minDist < pseudoSrcQ.top().dist))
		{
			Window win = winQ.top(); winQ.pop();
			if (win.level > mesh->m_nFace) continue;
			PropagateWindow(win);
		}
		else if (!pseudoSrcQ.empty() && (winQ.empty() || winQ.top().minDist >= pseudoSrcQ.top().dist))
		{
			PseudoWindow pseudoWin = pseudoSrcQ.top(); pseudoSrcQ.pop();
			if (pseudoWin.level >= mesh->m_nFace) continue;
			GenSubWinsForPseudoSrc(pseudoWin);
		}
	}
}

void ICH::BuildGeodesicPathTo(unsigned vertId)
{
	// TODO: build geodesic path from vertex vertId to source
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

			unsigned opVert = mesh->m_pEdge[mesh->m_pVertex[srcId].m_piEdge[j]].m_iVertex[1];
			vertInfos[opVert].birthTime = 0;
			vertInfos[opVert].dist = mesh->m_pEdge[mesh->m_pVertex[srcId].m_piEdge[j]].m_length;
			vertInfos[opVert].enterEdge = mesh->m_pEdge[mesh->m_pVertex[srcId].m_piEdge[j]].m_iTwinEdge;

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
			(directDist+win.pseudoSrcDist - splitInfos[e0].dist) / splitInfos[e0].dist > RELATIVE_ERROR)
		{
			hasLeftChild = splitInfos[e0].x < interX;
			hasRightChild = !hasLeftChild;
			/*cout << "Filter 1 works..." << endl;*/
		}
		else
		{
			splitInfos[e0].dist = directDist + win.pseudoSrcDist;
			splitInfos[e0].pseudoSrcId = win.pseudoSrcId;
			splitInfos[e0].srcId = win.srcID;
			splitInfos[e0].level = win.level;
			splitInfos[e0].x = l0 - interX;

			if (directDist + win.pseudoSrcDist < vertInfos[opVert].dist)
			{
				++vertInfos[opVert].birthTime;
				vertInfos[opVert].dist = directDist + win.pseudoSrcDist;
				vertInfos[opVert].enterEdge = e0;
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

	if (hasLeftChild) winQ.push(leftChildWin);
	if (hasRightChild) winQ.push(rightChildWin);

}

void ICH::GenSubWinsForPseudoSrc(const PseudoWindow &pseudoWin)
{
	unsigned startEdge, endEdge;
	if (mesh->m_pEdge[vertInfos[pseudoWin.vertID].enterEdge].m_iVertex[0] == pseudoWin.vertID)
		GenSubWinsForPseudoSrcFromPseudoSrc(pseudoWin, startEdge, endEdge);
	else if (mesh->m_pEdge[mesh->m_pEdge[vertInfos[pseudoWin.vertID].enterEdge].m_iNextEdge].m_iVertex[1] == pseudoWin.vertID)
		GenSubWinsForPseudoSrcFromWindow(pseudoWin, startEdge, endEdge);
	else assert(false);

	// generate windows
	while (startEdge != endEdge)
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

		startEdge = mesh->m_pEdge[mesh->m_pEdge[mesh->m_pEdge[startEdge].m_iNextEdge].m_iNextEdge].m_iTwinEdge;
	}

	// generate adjacent pseudo sources
	for (int i = 0; i < mesh->m_pVertex[pseudoWin.vertID].m_nValence; ++i)
	{
		unsigned adjEdge = mesh->m_pVertex[pseudoWin.vertID].m_piEdge[i];
		unsigned opVert = mesh->m_pEdge[adjEdge].m_iVertex[1];
		if (mesh->m_pAngles[opVert] < 2.0 * PI) continue;
		if (vertInfos[opVert].dist < pseudoWin.dist + mesh->m_pEdge[adjEdge].m_length) continue;

		vertInfos[opVert].dist = pseudoWin.dist + mesh->m_pEdge[adjEdge].m_length;
		++vertInfos[opVert].birthTime;
		vertInfos[opVert].enterEdge = mesh->m_pEdge[adjEdge].m_iTwinEdge;

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
	while (angle0 < PI || curEdge == -1)
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
	while (angle1 < PI || curEdge == -1)
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
	while (angle0 < PI || curEdge == -1)
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
	while (angle1 < PI || curEdge == -1)
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
		(win.pseudoSrcDist + (src2D - B).length()) / (vertInfos[v1].dist + win.b1) > RELATIVE_ERROR)
	{
		/*cout << "Filter 2 works..." << endl;*/
		return false;
	}
	if (win.pseudoSrcDist + (src2D - A).length() > vertInfos[v2].dist + l0 - win.b0 && 
		(win.pseudoSrcDist + (src2D - A).length()) / (vertInfos[v2].dist + l0 - win.b0) > RELATIVE_ERROR)
	{
		/*cout << "Filter 2 works..." << endl;*/
		return false;
	}
	if (isLeftChild)
	{
		if (win.pseudoSrcDist + (src2D - A).length() > vertInfos[v3].dist + (p3 - A).length() && 
			(win.pseudoSrcDist + (src2D - A).length()) / (vertInfos[v3].dist + (p3 - A).length()) > RELATIVE_ERROR)
		{
			/*cout << "Filter 2 works..." << endl;*/
			return false;
		}
	}
	else
	{
		if (win.pseudoSrcDist + (src2D - B).length() > vertInfos[v3].dist + (p3 - B).length() && 
			(win.pseudoSrcDist + (src2D - B).length()) / (vertInfos[v3].dist + (p3 - B).length()) > RELATIVE_ERROR)
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