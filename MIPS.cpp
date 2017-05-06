#include "MIPS.h"
#include <math.h>
#include <time.h>
#include <Eigen/SVD>

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;

#define pi 3.1415926
//Questions
//how to get my_mesh_? Answer : you don't need to get my_mesh_,you just write a method.
MIPS::MIPS()
{
}


MIPS::~MIPS()
{
	if (my_mesh != NULL)
	{
		delete my_mesh;
		my_mesh = NULL;
	}
	cout << "������������" << endl;
}

MIPS::MIPS(MyMesh *my_mesh_)
{
	my_mesh = my_mesh_;
}

double MIPS::sig_area_per_tri(const OpenMesh::FaceHandle &_fh)
{
	auto it = my_mesh->fv_begin(_fh);
	
	auto it1 = it;
	auto point1 = my_mesh->point(*it1);

	auto it2 = ++it;
	auto point2 = my_mesh->point(*it2);

	auto it3 = ++it;
	auto point3 = my_mesh->point(*it3);

	auto area = (point2 - point1) % (point3 - point2);
	//cout << area[0] << " " << area[1] << " " << area[2] << endl;
	if (area[2] > 0)
		return 1.;
	else
		return -1.;
}

double MIPS::sig_per_vert(const OpenMesh::VertexHandle &_vh)
{
	//double signal = 1.;
	for (auto it = my_mesh->vf_begin(_vh); it != my_mesh->vf_end(_vh); ++it)	
		if (sig_area_per_tri(*it) == -1.)
			return -1.;
	
	return 1.;
	//return signal;
}

void MIPS::the_inv_matrix()
{
	for (auto it = my_mesh->faces_begin(); it != my_mesh->faces_end(); it++)
	{
		//�������Լ�����������󣬶��ǻ���ͬһ���ֲ�����
		auto it1 = my_mesh->fv_begin(*it);
		auto n = my_mesh->normal(*it);//Ϊʲô��������ֱ�����أ���û��request��updateѽ

		auto point1 = my_mesh->point(*it1);
		++it1;
		auto point2 = my_mesh->point(*it1);
		++it1;
		auto point3 = my_mesh->point(*it1);

		auto x_axis = (point2 - point1).normalize();
		auto y_axis = n % x_axis;

		// x11 = 0. x12 = 0.
		// x22 = 0.
		auto x21 = (point2 - point1) | x_axis;
		auto x31 = (point3 - point1) | x_axis;
		auto x32 = (point3 - point1) | y_axis;

		Matrix2d m1;
		/*m1(0, 0) = -x21;
		m1(0, 1) = -x31;
		m1(1, 0) = 0.;
		m1(1, 1) = -x32;*/
		//Matrix2d m_inv = m1.inverse(); //���⣺ֱ��ʹ���治���ԣ���֪ԭ��
		//m1 = m1.inverse();
		//auto t =m.inverse();
		//Matrix2d mt = m.transpose();
		m1(0, 0) = -x32;
		m1(0, 1) = x31;
		m1(1, 0) = 0.;
		m1(1, 1) = -x21;
		m1 = m1 / (x21*x32);
		v_matrix1.push_back(m1);

		Matrix2d m2;
		m2(0, 0) = 0.;
		m2(0, 1) = -x21;
		m2(1, 0) = x32;
		m2(1, 1) = x21 - x31;
		m2 = m2 / (x21*x32);
		v_matrix2.push_back(m2);

		Matrix2d m3;
		m3(0, 0) = x32;
		m3(0, 1) = x21 - x31;
		m3(1, 0) = -x32;
		m3(1, 1) = x31;
		m3 = m3 / (x21*x32);
		v_matrix3.push_back(m3);

		//ÿһ�������̤ʵ
	}
	/*for (int i = 0; i < v_matrix1.size(); i++)
	{
	cout << v_matrix1[i](0, 0) << "  " << v_matrix1[i](0, 1) << "  " << v_matrix1[i](1, 0) << "  " << v_matrix1[i](1, 1) << endl;
	cout << v_matrix2[i](0, 0) << "  " << v_matrix2[i](0, 1) << "  " << v_matrix2[i](1, 0) << "  " << v_matrix2[i](1, 1) << endl;
	cout << v_matrix3[i](0, 0) << "  " << v_matrix3[i](0, 1) << "  " << v_matrix3[i](1, 0) << "  " << v_matrix3[i](1, 1) << endl;
	cout << " ******************************************* "<< endl;
	}*/
	/*cout << "the number of the vertex:" << my_mesh->n_faces() << endl;
	cout << "the size of v1 : " << v_matrix1.size() << endl;
	cout << "the size of v2 : " << v_matrix2.size() << endl;
	cout << "the size of v3 : " << v_matrix3.size() << endl;*/
}

//////////
//////////������1

void MIPS::uniform_para()
{
	the_inv_matrix();
	//������һ���Ŀ����Ϊ�˵�������ĺ������Բ������溯���Ĺ���
	//��һ����ʵ��ҲӦ�ý���

	////////project the boundary point to the unit circle

	//find the first boundary half_edge
	auto h_it = my_mesh->halfedges_begin();
	while (!my_mesh->is_boundary(*h_it))
		h_it++;

	//find the number of the boundary vertex(half_edge)
	int num = 0;
	auto h_start = *h_it;
	auto h_h = h_start;
	do
	{
		h_h = my_mesh->next_halfedge_handle(h_h);
		++num;
	} while (h_h != h_start);
	//num--;
	//num֮ǰȷʵ����������

	//project the boundary point to the xy_plane e.g set z = 0
	h_h = h_start;
	OpenMesh::VertexHandle hv;
	
	//��ֵ�����
	for (int i = 0; i < num; i++)
	{
		hv = my_mesh->to_vertex_handle(h_h);
		my_mesh->set_point(hv, OpenMesh::Vec3f(cos((- 2 * pi*i) / num), sin((- 2 * pi*i) / num), 0.));
		h_h = my_mesh->next_halfedge_handle(h_h);
	}

	//��ӡ���еĳ����������ж�
	/*h_h = h_start;
	do
	{
	std::cout << my_mesh->calc_edge_length(h_h) << std::endl;
	h_h = my_mesh->next_halfedge_handle(h_h);
	} while (h_h != h_start);*/

	//Ϊʲô����ĵ�ķֲ����Ǿ��ȵ�?
	//ԭ�򣺶������Ŀ������һ��,num����������


	//project the inner point to the xy_plane

	//fill A and bx,by
	std::vector<T> tripletList;
	OpenMesh::VectorT<double, 3> v_point;

	int m = my_mesh->n_vertices();
	SpMat A(m, m);
	VectorXd bx, by, x, y;
	//����ϸ�ڣ������֮ǰ��Xf�����ϸ����ĺ�ϸ������ѽ
	bx.resize(m);
	by.resize(m);

	//fill Triplet
	for (auto it1 = my_mesh->vertices_begin(); it1 != my_mesh->vertices_end(); ++it1)
	{
		auto it1_h = it1.handle();
		int id1 = it1_h.idx();
		bx[id1] = 0.;
		by[id1] = 0.;
		if (my_mesh->is_boundary(it1_h))
		{
			tripletList.push_back(T(id1, id1, 1));
			v_point = my_mesh->point(it1_h);
			bx[id1] = v_point[0];
			by[id1] = v_point[1];
		}
		else
		{
			tripletList.push_back(T(id1, id1, -(int)my_mesh->valence(it1_h)));
			/*	std::cout <<my_mesh->valence(it1_h)<<std::endl;
				std::cout << "hhhh"<<my_mesh->valence(it1_h) << std::endl;*/
			for (auto it2 = my_mesh->vv_begin(it1_h); it2 != my_mesh->vv_end(it1_h); ++it2)
			{
				auto it2_h = it2.handle();
				int id2 = it2_h.idx();
				tripletList.push_back(T(id1, id2, 1));
				//Ϊʲô�����һ����˵Ĵ���-1 ���� 1������Ļ��� -my_mesh->valence(it1_h)
				//����϶���tripletList�ķ�ʽ�й�ϵ
				//ԭ��my_mesh->valence()��������uint,��Ҫ��������ת��
			}
		}
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end());

	//solver Ax = b
	SparseLU<SparseMatrix<double> > solver;
	solver.compute(A);
	x = solver.solve(bx);
	y = solver.solve(by);

	//update all of the point
	for (auto it = my_mesh->vertices_begin(); it != my_mesh->vertices_end(); ++it)
	{
		int id = it.handle().idx();
		my_mesh->set_point(it.handle(), OpenMesh::Vec3f(x[id], y[id], 0));
	}

	my_mesh->request_face_normals();
	my_mesh->update_face_normals();
	for (auto it = my_mesh->faces_begin(); it != my_mesh->faces_end(); it++)
	{
		auto normal = my_mesh->normal(it.handle());
		//��������Щ����������Ķ���Щ�ĵ�
		//��δ�ӡ����normal��ֵ��������֮������normal��ֵӦ����һ���ģ�Ϊʲô��������
		//ԭ����û�ж�֮����request��update
		cout << "normal : " << normal[0] << " " << normal[1] << " " << normal[2] << " " << endl;
	}
}

double MIPS::detM(Matrix2d &m)
{
	double a = m(0, 0), b = m(0, 1), c = m(1, 0), d = m(1, 1);
	return a*d - b*c;
}

//����ÿ�������ζ�Ӧ��Jt����_vhΪδ֪��
Matrix2d MIPS::Jt_per_tri(const OpenMesh::FaceHandle &_fh, const OpenMesh::VertexHandle &_vh)
{
	auto it = my_mesh->fv_begin(_fh);
	//�������������أ����ڽ����
	auto it1 = it;
	auto point1 = my_mesh->point(*it1);

	auto it2 = ++it;//���ǰ++�ͺ�++�ǲ�һ����
	auto point2 = my_mesh->point(*it2);

	auto it3 = ++it;
	auto point3 = my_mesh->point(*it3);

	int id = _fh.idx();

	Matrix2d m;
	if (_vh == it1.handle())//���ǰ���֮ǰ����Ļ������Ӧ����_vh = it3;�����ᵼ����ȫ��Ӧ���󣬶�Ӧ����
	{
		m(0, 0) = (point1 - point2)[0];
		m(0, 1) = (point1 - point3)[0];
		m(1, 0) = (point1 - point2)[1];
		m(1, 1) = (point1 - point3)[1];

		//cout << "detm : " << detM(m) << " id : " << id << " inv_m : " << detM(v_matrix1.at(id)) << endl;
		m = m*v_matrix1.at(id);//��Ȼ�����idд���˵��id������һ����ɵ�Ƶ�����
	}
	else{
		if (_vh == it2.handle())
		{
			m(0, 0) = (point2 - point3)[0];
			m(0, 1) = (point2 - point1)[0];
			m(1, 0) = (point2 - point3)[1];
			m(1, 1) = (point2 - point1)[1];

			//cout << "detm : " << detM(m) << " id : " << id << " inv_m : " << detM(v_matrix2.at(id)) << endl;
			m = m*v_matrix2.at(id);		
		}
		else{
			m(0, 0) = (point3 - point1)[0];
			m(0, 1) = (point3 - point2)[0];
			m(1, 0) = (point3 - point1)[1];
			m(1, 1) = (point3 - point2)[1];
			//cout << "detm : " << detM(m) << " id : " << id << " inv_m : " << detM(v_matrix3.at(id)) << endl;
			m = m*v_matrix3.at(id);
		}
	}
	//cout << "m(0,0): " << m(0, 0) << " m(0,1): " << m(0, 1) << "m(1,0): " << m(1,0) << "m(1,1): " << m(1, 1) << endl;
	return m;
}

//Matrix2d MIPS::Jt_only_face(const OpenMesh::FaceHandle &)
//{
//
//}

//������Ľ�����Կ�������ط�����ر����أ�ԭ�����������ô��ˣ�����_fh��ƫƫʹ�õĵ����������

//////////
//////////������1

double MIPS::e_per_tri(const OpenMesh::FaceHandle &_fh, const OpenMesh::VertexHandle &_vh)
{
	Matrix2d m = Jt_per_tri(_fh, _vh);
	double a = m(0, 0), b = m(0, 1), c = m(1, 0), d = m(1, 1);
	double p = a*a + b*b + c*c + d*d;

	////�������Ҷ�djtȡ��һ������ֵ
	double djt = a*d - b*c;//����djtΪʲô�Ǹ��ģ�

	//cout << "djt : " << djt << endl;
	//cout << m(0, 0) << " " << m(0, 1) << " " << m(1, 0) << " " << m(1, 1);
	return p / djt;
}

//////////
//////////������1

Vector2d MIPS::grad_per_tri(const OpenMesh::FaceHandle &_fh, const OpenMesh::VertexHandle &_vh)
{
	//��Ӧ�������
	Matrix2d inv_m;
	auto it = my_mesh->fv_begin(_fh);
	auto it1 = it;
	auto it2 = ++it;
	auto it3 = ++it;

	if (_fh == it1.handle())
		inv_m = v_matrix1.at(_fh.idx());
	else
	{
		if (_fh == it2.handle())
			inv_m = v_matrix2.at(_fh.idx());
		else
			inv_m = v_matrix3.at(_fh.idx());
	}

	//��Ӧ��Jt
	Matrix2d m = Jt_per_tri(_fh, _vh);
	double a = m(0, 0), b = m(0, 1), c = m(1, 0), d = m(1, 1);
	double p = a*a + b*b + c*c + d*d;
	double djt = a*d - b*c;
	double djt2 = djt*djt;
	//double det_Jt = (m(0, 0)*m(1, 1) - m(0, 1)*m(1, 0)) * (m(0, 0)*m(1, 1) - m(0, 1)*m(1, 0));

	Vector2d v;
	v[0] = 0., v[1] = 0.;

	//Kx
	double Ka = 2 * a / djt - d*p / djt2;
	double ax = inv_m(0, 0) + inv_m(1, 0);
	v[0] += Ka*ax;

	double Kb = 2 * b / djt + c*p / djt2;
	double bx = inv_m(0, 1) + inv_m(1, 1);
	v[0] += Kb*bx;

	//Ky
	double Kc = 2 * c / djt + b*p / djt2;
	double cy = inv_m(0, 0) + inv_m(1, 0);
	v[1] += Kc*cy;

	double Kd = 2 * d / djt - a*p / djt2;
	double dy = inv_m(0, 1) + inv_m(1, 1);
	v[1] += Kd*dy;

	return v;
}

//////////
//////////������1

Matrix2d MIPS::inv_hessMat_per_tri(const OpenMesh::FaceHandle &_fh, const OpenMesh::VertexHandle &_vh)
{
	Matrix2d M;
	//m(0, 0) = 0., m(0, 1) = 0., m(1, 0) = 0., m(1, 1) = 0.; 

	//��Ӧ�������
	Matrix2d inv_m;
	auto it1 = my_mesh->fv_begin(_fh);
	auto it2 = it1++;
	auto it3 = it1++;
	if (_vh == it1.handle())//����֮ǰд���ǡ�_fh��ô��������������
		inv_m = v_matrix1.at(_fh.idx());
	else
	{
		if (_vh == it2.handle())
			inv_m = v_matrix2.at(_fh.idx());
		else
			inv_m = v_matrix3.at(_fh.idx());
	}

	//ax=cy , bx = dy
	double ax = inv_m(0, 0) + inv_m(1, 0);
	double cy = ax;
	double bx = inv_m(0, 1) + inv_m(1, 1);
	double dy = bx;

	//kaa,kab,kbb,kcc,kcd,kdd
	//��Ӧ��Jt
	Matrix2d m = Jt_per_tri(_fh, _vh);
	double a = m(0, 0), b = m(0, 1), c = m(1, 0), d = m(1, 1);
	double djt = a*d - b*c;

	double p = a*a + b*b + c*c + d*d;
	double djt2 = djt*djt;
	double djt3 = djt2*djt;

	double kaa = 2 / djt - 4 * a*d / djt2 + 2 * d*d*p / djt3;
	double kab = 2 * (a*c - b*d) / djt2 - 2 * c*d*p / djt3;
	double kac = 2 * (a*b - c*d) / djt2 - 2 * b*d*p / djt3;
	double kad = (-2 * d*d - 2 * a*a - p) / djt2 + 2 * a*d*p / djt3;

	double kba = kab;
	double kbb = 2 / djt + 4 * b*c / djt2 + 2 * c*c*p / djt3;
	double kbc = (2 * c*c + 2 * d*d + p) / djt2 + 2 * b*c*p / djt3;
	double kbd = 2 * (c*d - a*b) / djt2 - 2 * a*c*p / djt3;

	double kcc = 2 / djt + 4 * b*c / djt2 + 2 * b*b*p / djt3;
	double kcd = 2 * (b*d - a*c) / djt2 - 2 * a*b*p / djt3;

	double kdc = kcd;
	double kdd = 2 / djt - 4 * a*d / djt2 + 2 * a*a*p / djt3;

	//kxx,kxy,kyy
	//��Ȼ�������þ�����ʽ���������ֵ
	double kxx = (kaa * ax + kab * bx) * ax + (kba * ax + kbb * bx) * bx;
	double kxy = (kac * cy + kad * dy) * ax + (kbc * cy + kbd * dy) * bx;
	double kyx = kxy;
	double kyy = (kcc * cy + kcd * dy) * cy + (kdc * cy + kdd * dy) * dy;

	//��ʼ��M
	M(0, 0) = kxx, M(0, 1) = kxy, M(1, 0) = kyx, M(1, 1) = kyy;
	// 	double detM = kxx*kyy - kxy*kyx;
	// 	M(0, 0) = kyy, M(0, 1) = -kxy, M(1, 0) = -kyx, M(1, 1) = kxx;
	// 	M /= detM;
	return M;
}


//����ÿ������ľֲ�����
//Ϊ�ˣ��ȼ���ÿ��������������Ȼ����ӣ�
//��Ϊ�˼���ÿ�������ε�����������Ҫ����ÿ�������εı任����Jt
double MIPS::local_e(const OpenMesh::VertexHandle &_vh)
{
	double energy = 0.;
	//int num = 0;
	for (auto it = my_mesh->vf_begin(_vh); it != my_mesh->vf_end(_vh); ++it)
		energy += e_per_tri(it.handle(), _vh);

	return energy;
}

//////////
//////////������1

//compute gradient of local energy of each point
Vector2d MIPS::grad_of_le(const OpenMesh::VertexHandle &_vh)
{
	Vector2d v;
	v[0] = 0., v[1] = 0.;
	for (auto it = my_mesh->vf_begin(_vh); it != my_mesh->vf_end(_vh); it++)
		v += grad_per_tri(it.handle(), _vh);
	return v;
}

//////////
//////////������1

Matrix2d MIPS::inv_hessMat_le(const OpenMesh::VertexHandle &_vh)
{
	Matrix2d m;
	m(0, 0) = 0., m(0, 1) = 0., m(1, 0) = 0., m(1, 1) = 0.;
	for (auto it = my_mesh->vf_begin(_vh); it != my_mesh->vf_end(_vh); it++)
		m = m + inv_hessMat_per_tri(it.handle(), _vh);

	double a = m(0, 0), b = m(0, 1), c = m(1, 0), d = m(1, 1);
	m(0, 0) = d, m(0, 1) = -b, m(1, 0) = -c, m(1, 1) = a;

	m /= (a*d - b*c);
	return m;
}

double MIPS::maxi_step_per_tri(const OpenMesh::FaceHandle &_fh, const OpenMesh::VertexHandle &_vh, const Vector2d &d)
{
	auto it = my_mesh->fv_begin(_fh);
	auto it1 = it;
	auto it2 = ++it;
	auto it3 = ++it;

	auto it1_h = it1.handle();
	auto it2_h = it2.handle();
	auto it3_h = it3.handle();

	auto point1 = my_mesh->point(it1_h); 
	//cout << "point1[0]: " << point1[0] << "point1[1]: "<< point1[1] << endl;
	auto point2 = my_mesh->point(it2_h);
	//cout << "point2[0]: " << point2[0] << "point2[1]: " << point2[1] << endl;
	auto point3 = my_mesh->point(it3_h);
	//cout << "point3[0]: " << point3[0] << "point3[1]: " << point3[1] << endl;

	if (_vh == it1_h)
	{
		double fen_mu = (point2 - point3)[1] * d[0] + (point3 - point2)[0] * d[1];
		//cout << "fen_mu: " << fen_mu << endl;
		if (!fen_mu)
			return 0.;
		else
		{
			double fen_zi = (point1 - point3)[0] * (point1 - point2)[1] - (point1 - point3)[1] * (point1 - point2)[0];
			//cout << "fen_zi : " << fen_zi << endl;
			return fen_zi / fen_mu;
		}
	}
	else
	{
		if (_vh == it2_h)
		{
			double fen_mu = (point3 - point1)[1] * d[0] + (point1 - point3)[0] * d[1];
			//cout << "fen_mu: " << fen_mu << endl;
			if (!fen_mu)//����жϺ���������
				return 0.;
			else
			{
				double fen_zi = (point2 - point1)[0] * (point2 - point3)[1] - (point2 - point1)[1] * (point2 - point3)[0];
				//cout << "fen_zi : " << fen_zi << endl;
				return fen_zi / fen_mu;
			}
		}
		else
		{
			double fen_mu = (point1 - point2)[1] * d[0] + (point2 - point1)[0] * d[1];
			//cout << "fen_mu: " << fen_mu << endl;
			if (!fen_mu)
				return 0.;
			else
			{
				double fen_zi = (point3 - point2)[0] * (point3 - point1)[1] - (point3 - point2)[1] * (point3 - point1)[0];
				//cout << "fen_zi : " << fen_zi << endl;
				return fen_zi / fen_mu;
			}
		}
	}
}

//�о������������������⣬�������ĸ�

double MIPS::maxi_step_per_vert(const OpenMesh::VertexHandle &_vh, const Vector2d &d)
{
	////���õ�һ��������ʼ��step
	//auto it = my_mesh->vf_begin(_vh);
	//double step = maxi_step_per_tri(*it, _vh, d);
	//++it;
	//����㷨���ԣ�

	//����дһ������㷨
	if (my_mesh->is_boundary(_vh))
		return 0.;
	vector<double> v_step;
	for (auto it = my_mesh->vf_begin(_vh); it != my_mesh->vf_end(_vh); ++it)
	{
		double it_step = maxi_step_per_tri(*it, _vh, d);
		//cout << "it_step" << it_step << endl;
		//ԭ���п���ȫ����it_stepȫ��Ϊ��������϶�������ĳ�������
		if (it_step > 0.)
			v_step.push_back(it_step);
	}
	//cout << "v_step : " << v_step.size() << endl;
	if (v_step.size() == 0)//�������ǵ��Ǳ߽��ϵĵ������󲽳��п��ܶ��Ǹ��ģ�����Ҫ���ؿ���
		return 1.;
	else{
		double step = v_step[0];
		for (int i = 1; i < v_step.size(); i++)
		{
			if (step > v_step.at(i))
				step = v_step.at(i);
		}
		return step;
	}	
}

void MIPS::newton_Method()
{
	//�Ƚ��е���һ�Σ�Ȼ����˼��������Qt�Ի�����н�����ȷ�������Ĵ�����
	time_t t_start , t_end;
	t_start = time(NULL);
	uniform_para();
	for (int i = 0; i < 10000; ++i)
	{
		for (auto it = my_mesh->vertices_begin(); it != my_mesh->vertices_end(); ++it)
		{
			//cout << "signal : " << sig_per_vert(*it) << endl;
			/*if (sig_per_vert(*it) == -1.)
				cout << "this exist a face ,and a very beautiful face" << endl;*/
			const auto it_h = it.handle();//���ָ�ֵ����Ҫע��

			//����ÿ�������1����ľֲ�����;
			double k_original = local_e(it_h);
			//cout << "k_original :" << k_original << endl;
			//����ÿ��������ݶ�
			Vector2d g = grad_of_le(it_h);
			//_vg.push_back(g);
			//����ÿ�������Hessian�������
			Matrix2d inv_H = inv_hessMat_le(it_h);
			//  		_vh.push_back(inv_H);
			//  		cout << inv_H(0, 0) << " " << inv_H(0, 1) << " " << inv_H(1, 0) << " " << inv_H(1, 1) << endl;

			//  		_vd.push_back(d);
			

			//auto d = -g;//����100����ʱ24s?
			Vector2d d = -inv_H*g;//�������һ��Ҳ��24s?
			//double a = maxi_step_per_vert(*it, d);
			//��ʱd�ķ���Ϊ-g
			//cout << "a:= " << a << endl;
			//  		_va.push_back(a);
			//  		cout << "a:= " << a << endl;	


			//double a = 0.1;//��������ԣ�Ϊʲô�أ�
			double a = 0.1;//�����㷨��������õ�a��ֵΪ0.01��̫����
			//int i = 0;
			int id = it->idx();
			//cout << "id: " << id << endl;
			const auto ori_p = my_mesh->point(it_h);
			//handleû�䣬����point���ˡ�
			do{
				a /= 2;
				Vector2d ori_p_Eigen;
				ori_p_Eigen[0] = ori_p[0];
				ori_p_Eigen[1] = ori_p[1];
				Vector2d new_p;

				new_p = ori_p_Eigen + a*d;
				//
				my_mesh->set_point(*it, OpenMesh::Vec3f(new_p[0], new_p[1], 0));
				i++;
				//  	    cout << " i: " << i << std::endl;
				//  		cout << "k_original " << k_original << " " << "local_e(*it) : " << local_e(*it) << endl;
			} while (k_original < local_e(*it) || sig_per_vert(*it) == -1.);
		}
	}
	t_end = time(NULL);

	double delta_sum = 0.;
	for (auto it = my_mesh->faces_begin(); it != my_mesh->faces_end(); ++it)
	{
		//get the matrix
		auto _fh = it.handle();
		auto _vh = my_mesh->fv_begin(_fh).handle();
		Matrix2d m = Jt_per_tri(_fh, _vh);

		//acquire the SUV decomposition
		JacobiSVD<Matrix2d> svd(m);
		auto v = svd.singularValues();
		double sigmal_max = v[0];
		double sigmal_min = v[1];
		double delta = sigmal_max / sigmal_min;
		delta_val_per_tri.push_back(delta);
		delta_sum += delta;
		//cout << "delta : " << delta << endl;
	}
	double ave_delta = delta_sum / delta_val_per_tri.size();
	double maxi_mum = 0;
	double stan_error = 0.;
	for (int i = 0; i < delta_val_per_tri.size(); ++i)
	{
		double current = delta_val_per_tri[i];
		if (maxi_mum < current)
			maxi_mum = current;
		stan_error += (current - ave_delta)*(current - ave_delta);
		stan_error /= delta_val_per_tri.size();
		stan_error = sqrt(stan_error);
	}
	cout << "Average is : " << ave_delta << endl;
	cout << "Maximum is : " << maxi_mum << endl;
	cout << "Standard error is : " << stan_error << endl;

	printf("time: %.0f s\n", difftime(t_end, t_start));
}//�����е�2d�ĳ�2d