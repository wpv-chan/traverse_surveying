#pragma once
class traverse_survey
{
private:
	const double pi = 3.1415926;
public:
	traverse_survey();  //���߲����Ĺ��캯��
	~traverse_survey(); // ���߲�������������
	void h_d(double* d, double* m, double* s, int N, double* Decimal_degrees);//�ѹ۲�Ķȷ���ת���˻���ʹ��
	void transfrom_d_m_s(double Decimal_degrees); // ʮ���ƵĶ�ת��Ϊ�ȷ���
	void cal_trigo_function(double Decimal_degrees); // �������Ǻ���ֵ
	void position_angle(double x, double y);  //��λ��
	void start_end_position(double A_x, double A_y, double B_x, double B_y, double& start_position);//����ʼ��ĩλ�����귽λ�ǡ�
	void sum_distance(int N, double* distance, double& sum);//���ߵ��ܳ�
	void Angle_position(double start_positon, int N, double* Decimal_degrees, double* every_position);//��ÿһվ�����귽λ��
	void Correct_Angle_position(int N, double* every_position, double end_position, double* correct_every_position);//�������ÿһվ�����귽λ��
	void d_x_y(int N, double* distance, double* correct_every_position, double* d_x, double* d_y);//��d_x_y��ֵ
	void Cf_x(int N, double* d_x, double A_x, double C_x, double& f_x);//��f_x
	void Cf_y(int N, double* d_y, double A_y, double C_y, double& f_y); //��f_y
	void Cf_s(double sum, double f_x, double f_y); //��ȫ���պϲ��ȫ����Ապϲ�
	void Correct_d_x_y(int N, double* d_x, double* d_y, double* distance, double sum, double f_x, double f_y, double* correct_d_x, double* correct_d_y);//�����֮���d_x_y
	void x_y(int N, double A_x, double A_y, double* correct_d_x, double* correct_d_y, long double* lx, long double* ly);//��x,y;

};

