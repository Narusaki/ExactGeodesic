#include <Mesh.h>
#include <ICH.h>
#include <fstream>
#include <ctime>

using namespace std;

int main(int argc, char **argv)
{
	if (argc < 5)
	{
		cout << "USAGE: [.exe] [in.obj] [src0] [... [srcN]] [out.obj] [out.dist]" << endl;;
		return -1;
	}
	CMesh *mesh = new CMesh();
	if (!mesh->Load(argv[1]))
	{
		cout << "Cannot load mesh " << argv[1] << endl;
		return -2;
	}

	ICH *ich = new ICH();

	ich->AssignMesh(mesh);
	for (int i = 2; i < argc - 2; ++i)
		ich->AddSource(atoi(argv[i]));
	cout << "Executing ..." << endl;
	
	clock_t start = clock();
	ich->Execute();
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
		cout << "Cannot open output file " << argv[argc - 1] << endl;
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
	return 0;
}