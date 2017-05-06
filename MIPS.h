#pragma once
#include <iostream>
#include <OpenMesh\Core\Mesh\PolyConnectivity.hh>
#include <OpenMesh\Core\IO\MeshIO.hh>
#include <OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh>
#include <Eigen\SparseLU>
#include <vector>

using namespace std;
using namespace Eigen;

//using namespace OpenMesh;

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

class MIPS
{
public:
	MIPS();
	~MIPS();

public:
	MIPS(MyMesh *my_mesh_);

public:
	double sig_area_per_tri(const OpenMesh::FaceHandle &);
	double sig_per_vert(const OpenMesh::VertexHandle &);

	void the_inv_matrix();
	void uniform_para();
	double detM(Matrix2d &);

	Matrix2d Jt_per_tri(const OpenMesh::FaceHandle &, const OpenMesh::VertexHandle &);

	//����ֵ������ Matrix2d ����Ҫ��Ҫ��&�����أ�

	//�𰸣�
	
	//�������ֵ�����ڲ��������ͣ�������A GetA(void) ��дΪconst A & GetA(void)��ȷ�����Ч�ʡ�����ʱǧ��ǧ��ҪС�ģ�
	//һ��Ҫ����������������뷵��һ������ġ����������ǽ����ء��������Ϳ����ˣ������������
	double e_per_tri(const OpenMesh::FaceHandle &, const OpenMesh::VertexHandle &);
	Vector2d  grad_per_tri(const OpenMesh::FaceHandle &, const OpenMesh::VertexHandle &);
	Matrix2d inv_hessMat_per_tri(const OpenMesh::FaceHandle &, const OpenMesh::VertexHandle &);
	//FaceHandle ǰ��Ҫ��Ҫ����OpenMesh::�أ�
	//��Ҫ��,�����޶���


	double local_e(const OpenMesh::VertexHandle &);
	Vector2d grad_of_le(const OpenMesh::VertexHandle &);  
	double maxi_step_per_tri(const OpenMesh::FaceHandle &, const OpenMesh::VertexHandle &, const Vector2d &);
	double maxi_step_per_vert(const OpenMesh::VertexHandle &_vh, const Vector2d &);
	Matrix2d inv_hessMat_le(const OpenMesh::VertexHandle &);

public:
	void newton_Method();

private:
	MyMesh				  *my_mesh;
	//MyMesh              *my_mesh_again;
	vector<Matrix2d>	  v_matrix1;
	vector<Matrix2d>      v_matrix2;
	vector<Matrix2d>      v_matrix3;
	vector<Vector2d>	  _vv;
	vector<Vector2d>      _vg;
	vector<Vector2d>      _vd;
	vector<Matrix2d>      _vh;
	vector<double>         _va;

	vector<double>         delta_val_per_tri;
};

