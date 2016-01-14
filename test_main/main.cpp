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