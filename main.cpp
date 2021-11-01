//------------------------original main-----------------------------
//
//#define _USE_MATH_DEFINES
//
//#include "armadillo"
//
//#include <iostream>
//#include <string>
//#include <vector>
//#include <cmath>
//#include <tuple>
//#include <opencv2/opencv.hpp>
//
//
//using namespace std;
//using namespace arma;

#include "getCorresspondence.h"

class Point
{
private:
	double x, y, z;

public:

	Point(double nx=0, double ny =0,double nz=0):x(nx),y(ny),z(nz) {  }
	
	Point(const Point& rhs)
	{
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
	}
	Point& operator=(const Point& rhs)
	{
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;

		return *this;
	}
	~Point()
	{
		
	}
	const double getX() { return x; }
	const double getY() { return y; }
	const double getZ() { return z; }

	void setX(double num) { x = num; }
	void setY(double num) { y = num; }
	void setZ(double num) { z = num; }


	void setPoint(double nx, double ny, double nz)
	{
		setX(nx);
		setY(ny);
		setZ(nz);
	}
	
	const double getPoint() { }

	mat Point_to_mat(Point& rhs,string name)
	{
		
		
		
	}

	void ShowPoint()
	{
		cout << "x : " << this->x << " y : " << this->y << " z : " << this->z << endl;
	}
};
//
//struct St_correspondence_mat
//{
//	mat CadData;
//	mat measureed_data;
//	mat error_cal;
//	mat error;
//	mat Tr;
//	mat correspondence;
//	mat Translation;
//	mat Rotation;
//	mat err=eye(8,1);
//};

void fun_getLoadData(St_correspondence_mat *st,double row_num, double col_num, double row_distance, double col_distance)
{
	mat temp;	
	for (int i = 1; i <= row_num; i++) {
		for (int j = 1; j <= col_num; j++) {

			mat data(1,3);
			data.at(0) = ((i - 1) * row_distance);
			data.at(1) = ((j - 1) * col_distance);
			data.at(2) = 0;
			temp = join_cols(temp, data);			
		}
	}
	st->CadData = temp;


}


void point_to_mat(St_correspondence_mat *st,Point *rhs)
{
	mat temp;
	
	for (int j = 0; j < 8; j++) {
		for (int i = 0; i < 3; i++) {
			mat data(1, 3);
			data.at(i,j) = rhs->getX();
			data.at(i,j+1) = rhs->getY();
			data.at(i,j+2) = rhs->getZ();
			temp = join_cols(temp, data);
		}
	}

	st->measureed_data = temp;
}

void fun_getTransform(St_correspondence_mat *st, arma::mat M, arma::mat M0)
{
	//auto mean_M = arma::mean(XYZ);

	arma::mat XYZ = M;
	arma::mat Tr_XYZ = M0;
	arma::mat mean_M = arma::eye(1, 3);
	mean_M = arma::mean(XYZ);
	//cout << "mean_M : " << mean_M;
	auto numofPoints = XYZ.n_rows;
	st->Translation = mean_M - arma::mean(M0);
	//cout << "Translation : " << st->Translation << endl;

	mat XYZ_sub1 = XYZ.col(0);
	//XYZ_sub1.print("XYZ_sub1");

	auto Tr_XYZ_sub1 = mean_M.at(0);
	//cout << "Tr_XYZ_sub1 : " << Tr_XYZ_sub1 << endl;

	XYZ_sub1 = XYZ_sub1 - Tr_XYZ_sub1;
	//XYZ_sub1.print("1st col");

	auto XYZ_sub2 = XYZ.col(1);
	//XYZ_sub2.print("XYZ_sub2");

	auto Tr_XYZ_sub2 = mean_M.at(1);
	//cout << "Tr_XYZ_sub2 : " << Tr_XYZ_sub2 << endl;

	XYZ_sub2 = XYZ_sub2 - Tr_XYZ_sub2;
	//XYZ_sub2.print("2nd col");

	auto XYZ_sub3 = XYZ.col(2);
	//XYZ_sub3.print("XYZ_sub3");

	auto Tr_XYZ_sub3 = mean_M.at(2);
	//cout << "Tr_XYZ_sub2 : " << Tr_XYZ_sub3 << endl;

	XYZ_sub3 = XYZ_sub3 - Tr_XYZ_sub3;
	//XYZ_sub3.print("3rd col");

	auto res = arma::join_rows(XYZ_sub1, XYZ_sub2, XYZ_sub3);
	//res.print("res");

	auto M2 = res.t();
	//M2.print("M2");

	auto R = M2 * (Tr_XYZ) / numofPoints;
	//R.print("R");

	arma::mat U, V;
	arma::vec s;
	//arma::mat s = arma::eye(1, 3);
	arma::svd(U, s, V, R);
	

	/*U.print("U");
	s.print("s_vec");
	V.print("V");*/

	arma::mat S = arma::eye(3, 3);
	S.at(0) = s.at(0);
	S.at(4) = s.at(1);
	S.at(8) = s.at(2);
	//S.print("S_convert");

	if (arma::det(U*V) > 0) {
		st->Rotation = U * V.t();
		//st->Rotation.print("Rotation");

	}
	else {
		arma::mat bbc = arma::eye(3, 3);
		bbc.at(8) = -1;
		st->Rotation = U * bbc*V.t();
	}
	arma::mat TR = M.t() - st->Rotation * M0.t();
	//TR.print("TR");

	st->Translation = arma::mean(TR.t());
	//st->Translation.print("Translation");

	
}

mat findElement(double num, mat input)
{
	mat output(1, 2);
	for (int i = 0; i < input.n_rows; i++) {
		for (int j = 0; j < input.n_cols; j++) {
			if (input(i, j) == num)
			{
				output.at(0) = i;
				output.at(1) = j;
			}
		}
	}

	return output;
}


void fun_getCorrespondence(St_correspondence_mat *st,const mat o_CAD_data,const mat o_measured_data, int ci, int cj, int mi, int mj,const mat o_measured_V)
{
	auto virtual_pt_distance = 100;

	mat CAD_data = o_CAD_data;				//CAD_data.print("CAD_data");
	mat measured_data = o_measured_data;	//measured_data.print("measured_data");
	mat measuerd_v = o_measured_V;			//measuerd_v.print("measuerd_v");

	mat ith_CAD_data = CAD_data.row(ci);	//ith_CAD_data.print("ith_cad_data");
	mat jth_CAD_data = CAD_data.row(cj);	//jth_CAD_data.print("jth_cad_data");

	/*ith_CAD_data.print("ith cdata");
	jth_CAD_data.print("jth cdata");*/

	mat CAD_data_vecotr = (jth_CAD_data - ith_CAD_data) / arma::norm(jth_CAD_data - ith_CAD_data);
	//CAD_data_vecotr.print("cad_data_v");
	mat temp(1, 3);
	temp.fill(0);
	temp.at(2) = 1;
	//temp.print("temp");
	mat trd_CAD_data = ith_CAD_data + virtual_pt_distance * cross(temp,CAD_data_vecotr);
	mat M = join_cols(ith_CAD_data, jth_CAD_data, trd_CAD_data);
	//M.print("M");
	
	mat ith_MEA_data = measured_data.row(mi);
	//ith_MEA_data.print("ith_MEA_data");
	mat jth_MEA_data = measured_data.row(mj);
	//jth_MEA_data.print("jth_MEA_data");
	mat MEA_data_vector = (jth_MEA_data - ith_MEA_data) / norm(jth_MEA_data - ith_MEA_data);
	//MEA_data_vector.print("mea_Data_vector");
	mat trd_MEA_data = ith_MEA_data + virtual_pt_distance * cross(measuerd_v, MEA_data_vector);
	mat M0 = join_cols(ith_MEA_data, jth_MEA_data, trd_MEA_data);
	//M0.print("M0");
	
	fun_getTransform(st,M0, M);
	//st->Rotation.print("Rotation");
	//st->Translation.print("Translation");

	mat inv_Tr,ones;
	ones = eye(1, 4);
	ones = "0 0 0 1";
	//ones.print("ones");
	inv_Tr = join_rows(st->Rotation, st->Translation.t()); //, ones);
	inv_Tr = join_cols(inv_Tr, ones);
	//inv_Tr.print("inv_Tr");
	mat tRecal2,tRecal,ttemp;
	
	ones.reshape(1, CAD_data.n_rows);
	ones.fill(1);
	//ones.print("ones");

	ttemp = join_cols(CAD_data.t(), ones); //tRecal.print("tt");
	tRecal2 = inv_Tr * ttemp;
	tRecal = tRecal2.submat(span(0, 2), span(0, 7)).t();
	//tRecal.print("XYZ_recal");
	
	mat dist_affinity_matrix(8,8);
	dist_affinity_matrix.fill(NULL);
	//dist_affinity_matrix.print("dist_affinity_matrix");
	
	for (int i = 0; i < tRecal.n_rows; i++) {
		for (int j = 0; j < measured_data.n_rows; j++) {
			mat CAD_data_A = tRecal.row(i);
			mat MEA_data_A = measured_data.row(j);
			auto distance = sqrt(accu ( pow( ( CAD_data_A - MEA_data_A ), 2) ));
			dist_affinity_matrix(i, j) = distance;			
		}
	}
	//dist_affinity_matrix.print("dist_affinity_matrix");
	//dist_affinity_matrix.print("dist_affinity_matrix");
	auto max_Val = max(max(dist_affinity_matrix));
	mat copy_dist_affinity_matrix = dist_affinity_matrix;
	st->correspondence.reset();
	
	
	for (int k = 0; k < tRecal.n_rows * measured_data.n_rows; k++) {
		auto min_Val = arma::min(arma::min(copy_dist_affinity_matrix));
		if (min_Val < 10)
		{
			mat crtemp(1,2);
			
			//copy_dist_affinity_matrix.print("copy mat");
			crtemp = findElement(min_Val, copy_dist_affinity_matrix);
			//crtemp.print("find element");
			auto idx_Cad = crtemp.at(0);
			auto idx_Mea = crtemp.at(1);
			//cout << "row , col : " << idx_Cad << " , " << idx_Mea << endl;
			
			st->correspondence = join_cols(st->correspondence, crtemp);
			//st->correspondence.print("correspondence");

			mat cpy_row = copy_dist_affinity_matrix.row(idx_Cad);
			mat cpy_col= copy_dist_affinity_matrix.col(idx_Mea);
			cpy_row.fill(max_Val * 2);
			cpy_col.fill(max_Val * 2);
			
			copy_dist_affinity_matrix.row(idx_Cad) = cpy_row;
			copy_dist_affinity_matrix.col(idx_Mea) = cpy_col;			
			//copy_dist_affinity_matrix.print("after copy mat");



		}
		else
		{
			//cout << "break" << endl;
			break;
		}
	}

	if (st->correspondence.n_rows < 7) 
	{
		st->err.fill(100000);
		st->Tr = eye(4, 4);

		return;
	}

	mat M_CAD=eye(CAD_data.n_rows,CAD_data.n_cols);
	mat	M_MEA=eye(measured_data.n_rows,measured_data.n_cols);
	for (int m = 0; m < st->correspondence.n_rows; m++) {
		//CAD_data.print("cad_data");
		//measured_data.print("measured_data");
		//st->correspondence.print("corresspondence");
		mat c_row(1,3), m_row(1,3);
		

		M_CAD.row(m) = CAD_data.row(st->correspondence(m, 0)); //cout << CAD_data.row(st->correspondence(m, 0)) << endl;
		M_MEA.row(m) = measured_data.row(st->correspondence(m, 1));  //cout << measured_data.row(st->correspondence(m, 1)) << endl;

		/*for (int n = 0; n < 3; n++) {
			c_row.reset(); m_row.reset();
			c_row.at(m) = st->correspondence(m, 0);
			m_row.at(m) = st->correspondence(m, 1);			
		}*/
		
		//M_CAD = join_cols(M_CAD, c_row);
		//M_CAD.print("m_cad");
		//M_MEA = join_cols(M_MEA, m_row);
		//M_MEA.print("m_mea");
		//CAD_data.col(st->correspondence(m, 1)).print("test");


		//M_CAD.col(m) = CAD_data.col(st->correspondence(m, 1));
		//M_MEA.col(m) = CAD_data.col(st->correspondence(m, 2));


	}
	

	//M_CAD.print("m_cad"); M_MEA.print("m_mea");
	//st->Rotation.print("before rotation"); st->Translation.print("before translation");
	fun_getTransform(st, M_MEA, M_CAD);
	inv_Tr = join_rows(st->Rotation, st->Translation.t()); //, ones);
	//st->Rotation.print("after rotation"); st->Translation.print("after translation");
	ones.reshape(1, 4);
	ones = "0 0 0 1";
	inv_Tr = join_cols(inv_Tr, ones);

	//inv_Tr.print("inv_Tr");

	if (st->Rotation(0, 0) < 0)
	{
		st->err.fill(10000);
		st->Tr = eye(4, 4);

		return;
	}
	mat cal3_temp;
	ones.reshape(1, M_CAD.n_rows);
	ones.fill(1);
	cal3_temp = join_cols( M_CAD.t(),ones);

	mat xyz_re_cal3 = inv_Tr * cal3_temp;
	mat xyz_re_cal4 = xyz_re_cal3.submat(span(0, 2), span(0, 7)).t();

	//xyz_re_cal4.print("xyz_re_cal4");

	//M_MEA.print("m_mea"); xyz_re_cal4.print("xyz_re_cal4");

	mat dx(8, 1), dy(8, 1), dz(8, 1);
	dx = M_MEA.col(0) - xyz_re_cal4.col(0);		        //dx.print("dx");
	dy = M_MEA.col(1) - xyz_re_cal4.col(1);				//dy.print("dy");
	dz = M_MEA.col(2) - xyz_re_cal4.col(2);				//dz.print("dz");

	st->err =sqrt( pow(dx, 2) + pow(dy,2) + pow(dz,2));
	st->Tr = inv_Tr;
	
	
	
}

int main(void)
{
	
	////St_correspondence_mat *st;
	//getCorresspondence crs;
	//St_correspondence_mat *st = new St_correspondence_mat;
	///*getCorresspondence::St_correspondence_mat *st =
	// new getCorresspondence::St_correspondence_mat;
	//*/
	//
	////Point arr_p[8];
	////Point *p1 = new Point[8];
	////bool pointCheck = false;
	//int count = 1;
	//int i = 0; //count class	

	///*while (!pointCheck)
	//{
	//	if (count == 8) { pointCheck = true; }
	//	else {
	//		double mx = 0; double my = 0; double mz = 0;
	//		cin >> mx >> my >> mz;
	//		p1[i].setPoint(mx, my, mz);
	//		p1[i].ShowPoint();

	//		count++;
	//	}
	//}*/

	///*p1[0].setPoint(87.4638060000000 , - 107.085831000000,563);
	//p1[1].setPoint(8.38989300000000 ,- 107.284851000000	,568);
	//p1[2].setPoint(6.64999900000000 ,- 46.1396680000000	,570.000061000000);
	//p1[3].setPoint(90.9495700000000 ,- 45.4251480000000	,562.000061000000);
	//p1[4].setPoint(89.8627010000000	,10.1085700000000	,566.000061000000);
	//p1[5].setPoint(6.65303600000000	,14.6639440000000	,570.000061000000);
	//p1[6].setPoint(89.7951580000000	,70.3689040000000	,566.000061000000);
	//p1[7].setPoint(6.69246200000000	,73.2542340000000	,574);*/
	//
	//
	////st->measureed_data.load("measured_data.mat", raw_binary);
	////st->measureed_data.reshape(8, 3);

	////st->measureed_data = "87.4638060000000, -107.085831000000,	563;	8.38989300000000, -107.284851000000,	568;		6.64999900000000, -46.1396680000000,	570.000061000000;		90.9495700000000, -45.4251480000000,	562.000061000000;		89.8627010000000,	10.1085700000000,	566.000061000000;		6.65303600000000,	14.6639440000000,	570.000061000000;		89.7951580000000,	70.3689040000000,	566.000061000000;		6.69246200000000,	73.2542340000000,	574;";
	//st->measureed_data.load("measured_data.csv");
	////st->measureed_data.print("measured_data");
	////mat fromMatlab,fromCSV;
	////fromMatlab.load("measured_data.mat");
	////fromCSV.load("measured_data.csv");
	////fromMatlab.print("load .mat");
	////fromCSV.print("measured_data.csv");

	////fun_getLoadData(st,2,4,83,61);
	////fun_getLoadData(st,2, 4, 83, 61);
	//crs.fun_getLoadData(st, 2, 4, 83, 61);
	////st->CadData.print("CadData");
	//
	//mat U, V;
	//vec s;		
	//mat ones=eye(8,1);
	//ones.fill(1);

	//mat svd_input = join_rows(st->measureed_data, ones);
	////svd_input.print("svd_input");

	//
	//svd(U, s, V, svd_input);
	////U.print("U"); s.print("s"); V.print("V");
	//
	//

	//mat v1 = V.submat(0, 3, 2, 3).t(); //.t().print("v(3)");
	////v1.print("v1");
	//
	//auto v2 = arma::norm(v1);
	//mat measured_v(1, 3);
	//measured_v = v1 / v2;
	////measured_v.print("measured_v");

	//int cSize = st->CadData.n_rows;
	//int mSize = st->measureed_data.n_rows;
	//
	//mat err_temp(1,5);
	//st->error;

	//for (int cad_i = 0; cad_i < cSize; cad_i++) {
	//	for (int cad_j = 0; cad_j < cSize; cad_j++) {
	//		for (int mea_i = 0; mea_i < mSize; mea_i++) {
	//			for (int mea_j = 0; mea_j < mSize; mea_j++) {
	//				if (((cad_i == cad_j) || (mea_i == mea_j)))
	//				{
	//					//fun_getCorrespondence(st, st->CadData, st->measureed_data, cad_i, cad_j, mea_i, mea_j, measured_v);
	//				}
	//				else 
	//				{
	//					fun_getCorrespondence(st,st->CadData, st->measureed_data, cad_i, cad_j, mea_i, mea_j, measured_v);
	//					//error_cal = [ error_cal; [ CAD_i, CAD_j, MEA_i, MEA_j, mean(error) ] ];
	//					//auto meanErr = mean(st->err);
	//					mat ermean(8, 1);
	//					ermean = arma::mean(st->err);
	//					//ermean.print("ermean");
	//					//st->err.print("err");
	//					err_temp.at(0) = (cad_i);
	//					err_temp.at(1) = (cad_j);
	//					err_temp.at(2) = (mea_i);
	//					err_temp.at(3) = (mea_j);
	//					err_temp.at(4) = ermean(0);
	//					st->error_cal = join_cols(st->error_cal, err_temp);
	//					//st->error_cal.print("err_cal");
	//					
	//				}
	//			}
	//		}
	//	}
	//}
	////st->error_cal.print("err_cal");
	//mat index_min_error_cal;
	//mat test;
	//auto ans = arma::min(arma::min(st->error_cal.col(4)));
	//auto ans2 = arma::find(st->error_cal.col(4) == ans);
	////ans2.print("ans2");

	//index_min_error_cal = findElement(arma::min(st->error_cal.col(4)), st->error_cal.col(4));
	//int CAD_i = st->error_cal(index_min_error_cal(0), 0);
	//int CAD_j = st->error_cal(index_min_error_cal(0), 1);
	//int MEA_i = st->error_cal(index_min_error_cal(0), 2);
	//int MEA_j = st->error_cal(index_min_error_cal(0), 3);
	//fun_getCorrespondence(st, st->CadData, st->measureed_data, CAD_i, CAD_j, MEA_i, MEA_j, measured_v);

	//st->err.print("err");
	//st->Tr.print("Tr");
	//st->correspondence.print("correspondence");
	//mat total_err_mean(8, 1);
	//total_err_mean = mean(st->err);
	//total_err_mean.print("total err mean");
	////st->error_cal.col.print("err_cal");
	////cout << findElement(arma::min(st->error_cal.col(4)),st->error_cal.col(4)) <<endl;

	////cout << arma::min(st->error_cal.col(4)) << endl;


	//cout << endl << "end of process" << endl;
	getCorresspondence crs;
	crs.ShowMat();

	//delete st;
	//delete[] p1;

	return 0;
}