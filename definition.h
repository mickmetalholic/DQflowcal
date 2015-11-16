#ifndef _DEFINITION_H
#define	_DEFINITION_H

////////////////////////
///////�ṹ�嶨��///////
////////////////////////
/*****֧·*****/
struct Branch_Info
{
	int i,j;							//���ӽڵ�i��j
	double R,X,YK;						//֧·���衢�翹���ӵص���
};
/*****�����*****/
struct Generator_Info
{
	int i,j;							//���ӽڵ�i���ڵ�����j
	double P,Q,V;						//�ڵ��й����޹�������PV�ڵ�Ϊ�޹����ޣ�����ѹ������PV�ڵ�Ϊ��Ҫά�ֵ�ѹ ��Ϊ����
};
/*****����*****/
struct Load_Info
{
	int i,j;							//���ӽڵ�i���ڵ�����j
	double P,Q,V;						//�ڵ��й����޹�����ѹ
};
/*****�ڵ�*****/
struct Node_Info
{
	int num;							//�ڵ��
	int outnum;							//������
	int connect[20];					//��ýڵ������Ľڵ��
};
/*****PV�ڵ�*****/
struct PVNode_Info
{
	double V;							//������ѹ
	int i;								//�ڵ��
};
/*****�Խ�Ԫ��*****/
struct Yii_Info
{
	double G,B;                         //GΪʵ�� BΪ�鲿
};
/*****�ǶԽ�Ԫ��*****/
struct Yij_Info
{
	double G,B;                         //GΪʵ�� BΪ�鲿
	int j;								//�ǶԽ�Ԫ���к�
};
/*****���ӱ������Ǿ���*****/
struct U_Info
{
	double value;
	int j;
};
/*****�ڵ��ѹ*****/
struct NodalVoltage_Info
{
	double V,theta;						//�ڵ��ѹ��ֵV �Ƕ�theta
};
/*****�ڵ㹦��*****/
struct NodalPower_Info
{
	double P,Q;							//�ڵ��й�P �޹�Q
};


////////////////////////
////////��������////////
////////////////////////
void Read(char *filename);																	//���ݶ�ȡ
void NodeOpt();																				//�ڵ����Ż�
void YFormation();																			//���ɾ����γ�
void BFormation(Yii_Info *Yii,Yij_Info *Yij,int flag,int *NUsum,U_Info *U,double *D);		//�γ�B' B''�����ӱ�
void Initial(double V0,double theta0);														//���ڵ��ѹ��ֵ
void PQCal(int flag);																		//�ڵ㹦�ʼ���
void PQErrorCal(double *DI,int flag);																//�����������
void SolvEqu(U_Info *U,double *D,double *DI,int *NUsum);									//�����Է�����
void OutData();






#endif