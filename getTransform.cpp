// getTransform.cpp : 이 파일에는 'main' 함수가 포함됩니다. 거기서 프로그램 실행이 시작되고 종료됩니다.
//

#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <tuple>
#include <opencv2/opencv.hpp>
#include "armadillo"
using namespace std;
cv::Mat arma_to_cv(const arma::mat input)
{
	arma::mat in = input;
	cv::Mat out(in.n_cols, in.n_rows, CV_64F, in.memptr());

	return out.t();

}

arma::mat cv_to_arma(const cv::Mat input)
{
	arma::mat output = arma::eye(input.rows, input.cols);
	for (int i = 0; i < input.cols; i++)
		for (int j = 0; j < input.rows; j++)
			output.at(j, i) = input.at<double>(j, i);
	return output;
}

//struct St_rodrigues {
//	arma::mat R;
//	arma::mat out;
//	arma::mat omega;
//	arma::mat theta2;
//	arma::mat sinthetatheta;
//	arma::mat onecosthetatheta2;
//
//};

struct St_getTransfromedPoint {
	arma::mat Tr;
	arma::mat XYZ2A;
	arma::mat Tr_XYZ;
	arma::mat R;
	arma::mat out;
	arma::mat omega;
	arma::mat theta2;
	arma::mat sinthetatheta;
	arma::mat onecosthetatheta2;
};


struct St_getTransform {
	arma::mat Rotation;
	arma::mat Translation;
};

struct St_main {
	arma::mat U = "0,0,0,1";
};

arma::mat skew3(const arma::mat input)
{
	arma::mat s = arma::eye(3, 3);
	s.at(0, 0) = 0;
	s.at(0, 1) = -input.at(2);
	s.at(0, 2) = input.at(1);
	s.at(1, 0) = input.at(2);
	s.at(1, 1) = 0;
	s.at(1, 2) = -input.at(0);
	s.at(2, 0) = -input.at(1);
	s.at(2, 1) = input.at(0);
	s.at(2, 2) = 0;

	arma::mat dSdT = "0 0 0;0 0 1;0 -1 0 ;0 0 -1;0 0 0;1 0 0 ;0 1 0;-1 0 0; 0 0 0";

	//return tuple(s, dSdT);
	return s;
}

arma::mat rodrigues(St_getTransfromedPoint *st,const arma::mat in)
{
	//arma::mat inver;
	//input = arma::eye(3, 3);
	arma::mat input = in;
	//St_rodrigues  *st = new St_rodrigues;
	auto m = input.n_rows;
	auto n = input.n_cols;
	auto eps = 2.2204e-16;
	auto bigeps = (10e+5) * 2.2204e-16;
	bool bSuccess = false;

	if ((m == 1 && n == 3) || (m == 3 && n == 1))
	{
		auto theta = arma::norm(input);
		if (n == (input.n_rows)) {
			input = arma::inv(input);
		}
		st->omega = input;

		if (theta < sqrt(eps)*1e2) {
			st->theta2 = input * arma::pinv(input);
			auto thetatest = input * arma::inv(input);
			st->sinthetatheta = 1 - st->theta2 / 6;
			st->onecosthetatheta2 = 1 - st->theta2 / 24;
		}
		else {
			st->theta2 = theta * theta;
			st->onecosthetatheta2 = (1 - cos(theta)) / st->theta2;
			st->sinthetatheta = sin(theta) / theta;
		}

		auto omegav = skew3(st->omega);
		/*omegav.print("omegav");
		st.sinthetatheta.print("sinthetatheta");
		st.onecosthetatheta2.print("onecosthetatheta2");
		auto ans = arma::eye(3, 3);
		auto t = ans + omegav*st.sinthetatheta.at(0);
		t.print("ans");*/
		st->R = arma::eye(3,3) + omegav * st->sinthetatheta.at(0) + omegav * omegav * st->onecosthetatheta2.at(0);
		//st->R.print("R");
		st->out = st->R;
		
		auto R_ = st->R.t();
		//R_.print("R_");
	}

	
	return st->out;
}


void fun_getTransformedPoints(St_getTransfromedPoint *st_gform ,arma::mat v, double noiselevel, arma::mat input)
{
	
	double pi = M_PI;
	auto in = (v.at(0)*pi / 180, v.at(1)*pi / 180, v.at(2)*pi / 180);

	arma::mat v1 = arma::eye(1, 3);
	v1.at(0) = v.at(0)*pi / 180;
	v1.at(1) = v.at(1)*pi / 180;
	v1.at(2) = v.at(2)*pi / 180;

	auto v1_ = rodrigues(st_gform,v1);
	//v1_.print("v1_");

	arma::mat v2 = arma::eye(3, 1);
	v2.at(0) = v.at(3);
	v2.at(1) = v.at(4);
	v2.at(2) = v.at(5);
	//v2.print("v2");

	arma::mat rMat = arma::join_rows(v1_, v2);
	//rMat.print("rMat");

	arma::mat u = arma::eye(1, 4);
	u = "0, 0, 0, 1";
	//cout << "u : " << u << endl;

	st_gform->Tr = arma::join_cols(rMat, u);
	//st_gform->Tr.print("Tr");

	arma::mat xyz2a;
	arma::mat ones = arma::eye(1,input.n_rows);
	//input.t().print("xyz'");
	ones.fill(1);
	//ones.print("ones");
	arma::mat temp = arma::join_cols(input.t(), ones);
	//temp.print("temp");
	xyz2a = st_gform->Tr * temp;
	//xyz2a.print("xyz2a");
	

	arma::mat xyz2b;
	xyz2b = xyz2a.t();
	xyz2b = xyz2b.submat(0, 0, 4, 2);
	//xyz2b.print("xyz2b");

	arma::mat smat = arma::randn<arma::mat>(5, 3);
	auto noise = noiselevel * (smat / 3);
	//noise.print("noise");

	st_gform->Tr_XYZ = xyz2b + noise;
	//st_gform->Tr_XYZ.print("Tr_XYZ");

	//return input;
}

arma::mat fun_getTransform(St_getTransform *st,arma::mat M, arma::mat M0)
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
	
	auto XYZ_sub1 = XYZ.submat(0, 0, 4, 0);
	//XYZ_sub1.print("XYZ_sub1");

	auto Tr_XYZ_sub1 = mean_M.at(0);
	//cout << "Tr_XYZ_sub1 : " << Tr_XYZ_sub1 << endl;

	XYZ_sub1 = XYZ_sub1 - Tr_XYZ_sub1;
	//XYZ_sub1.print("1st col");

	auto XYZ_sub2 = XYZ.submat(0, 1, 4, 1);
	//XYZ_sub2.print("XYZ_sub2");

	auto Tr_XYZ_sub2 = mean_M.at(1);
	//cout << "Tr_XYZ_sub2 : " << Tr_XYZ_sub2 << endl;

	XYZ_sub2 = XYZ_sub2 - Tr_XYZ_sub2;
	//XYZ_sub2.print("2nd col");

	auto XYZ_sub3 = XYZ.submat(0, 2, 4, 2);
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

	arma::mat U ,V;
	arma::vec s;
	//arma::mat s = arma::eye(1, 3);
	arma::svd(U, s, V, R);
	
	//U.print("U");
	//s.print("s_vec");
	//V.print("V");
	
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
	
	return XYZ;	
}




int main()
{
	St_getTransfromedPoint *st_gTform = new St_getTransfromedPoint;
	St_getTransform *st_Trs = new St_getTransform;
	//double mem[] = { 100,200,100, 110,200,120, 80,220,200, 90,170,110,110,180,130 };
	arma::mat xyz = "100,200,100;110,200,120;80,220,200;90,170,110;110,180,130";
	//arma::mat xyz = arma::mat(mem, 5, 3);

	arma::mat diff = arma::eye(1,6);
	diff = "10 20 2 100 200 100";	
	
	double noise_level = 0.0;
	double pi = M_PI;
	std::cout << "xyz : " << std::endl << xyz << std::endl;
	std::cout << "diff : " << std::endl << diff << std::endl;
	std::cout << "noise level : " << std::endl;
	std::printf("%.1lf\n", noise_level);
		
	
	fun_getTransformedPoints(st_gTform,diff,noise_level,xyz); //return Tr,Tr_XYZ
	st_gTform->Tr.print("Tr");
	st_gTform->Tr_XYZ.print("Tr_xyz");

	auto cTr_XYZ = st_gTform->Tr_XYZ; //copy of Tr_XYZ
	cTr_XYZ.print("copy of Tr Xyz");
	fun_getTransform(st_Trs, xyz, cTr_XYZ);	 //return Rotation , Translation

	st_Trs->Rotation.print("Rotation");
	st_Trs->Translation.print("Translation");

	//inv_Tr = [ [ Rotation, Translation']; [ 0 0 0 1 ] ];
	arma::mat inv_Tr_temp = arma::join_rows(st_Trs->Rotation, st_Trs->Translation.t());
	inv_Tr_temp.print("inv_Tr_temp");
	St_main *mst = new St_main;
	

	//mst->U : 0001
	arma::mat inv_Tr=arma::join_cols(inv_Tr_temp, mst->U);
	inv_Tr.print("inv_Tr");

	arma::mat Tr_Cal = arma::inv(inv_Tr);
	Tr_Cal.print("Tr_Cal");

	arma::mat temp;
	arma::mat ones = arma::eye(1, 5);
	ones.fill(1);
	temp = arma::join_cols(cTr_XYZ.t(), ones);
	
	arma::mat XYZ_re_cal2 = inv_Tr * temp;
	XYZ_re_cal2.print("xyz_re_cal2");

	arma::mat XYZ_re_cal = XYZ_re_cal2.submat(0,0,2,4).t();
	XYZ_re_cal.print("xyz_re_cal");
	xyz.print("xyz");
	
	arma::mat err = arma::sqrt(
		arma::pow((xyz.col(0) - XYZ_re_cal.col(0)), 2) //  (x-x').^2
	  + arma::pow((xyz.col(1) - XYZ_re_cal.col(1)), 2) //+ (y-y').^2
      + arma::pow((xyz.col(2) - XYZ_re_cal.col(2)), 2) //+ (z-z').^2
	);

	std::cout << "xyz.col(0) - XYZ_re_cal.col(0) : " << std::endl << pow(xyz.col(0) - XYZ_re_cal.col(0) ,2) << endl;
	std::cout << "xyz.col(1) - XYZ_re_cal.col(1) : " << std::endl << pow(xyz.col(1) - XYZ_re_cal.col(1) ,2) << endl;
	std::cout << "xyz.col(2) - XYZ_re_cal.col(2) : " << std::endl << pow(xyz.col(2) - XYZ_re_cal.col(2) ,2)<< endl;
	
	err.print("err");

	/*arma::mat e;
	e.randn(4,4);
	e.print("e");
	e.i().print("e'");
	cout << e * e.i() << endl;*/
	
	
	
	

	
	delete st_gTform;
	delete st_Trs;
	delete mst;
	
	
	return 0;
	
}

// 프로그램 실행: <Ctrl+F5> 또는 [디버그] > [디버깅하지 않고 시작] 메뉴
// 프로그램 디버그: <F5> 키 또는 [디버그] > [디버깅 시작] 메뉴

// 시작을 위한 팁: 
//   1. [솔루션 탐색기] 창을 사용하여 파일을 추가/관리합니다.
//   2. [팀 탐색기] 창을 사용하여 소스 제어에 연결합니다.
//   3. [출력] 창을 사용하여 빌드 출력 및 기타 메시지를 확인합니다.
//   4. [오류 목록] 창을 사용하여 오류를 봅니다.
//   5. [프로젝트] > [새 항목 추가]로 이동하여 새 코드 파일을 만들거나, [프로젝트] > [기존 항목 추가]로 이동하여 기존 코드 파일을 프로젝트에 추가합니다.
//   6. 나중에 이 프로젝트를 다시 열려면 [파일] > [열기] > [프로젝트]로 이동하고 .sln 파일을 선택합니다.
