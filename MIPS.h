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

	//返回值类型中 Matrix2d 后面要不要加&符号呢？

	//答案：
	
	//如果返回值不是内部数据类型，将函数A GetA(void) 改写为const A & GetA(void)的确能提高效率。但此时千万千万要小心，
	//一定要搞清楚函数究竟是想返回一个对象的“拷贝”还是仅返回“别名”就可以了，否则程序会出错。
	double e_per_tri(const OpenMesh::FaceHandle &, const OpenMesh::VertexHandle &);
	Vector2d  grad_per_tri(const OpenMesh::FaceHandle &, const OpenMesh::VertexHandle &);
	Matrix2d inv_hessMat_per_tri(const OpenMesh::FaceHandle &, const OpenMesh::VertexHandle &);
	//FaceHandle 前面要不要加上OpenMesh::呢？
	//需要加,这是限定符


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

