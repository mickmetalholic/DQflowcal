#ifndef _DEFINITION_H
#define	_DEFINITION_H

////////////////////////
///////结构体定义///////
////////////////////////
/*****支路*****/
struct Branch_Info
{
	int i,j;							//连接节点i、j
	double R,X,YK;						//支路电阻、电抗、接地导纳
};
/*****发电机*****/
struct Generator_Info
{
	int i,j;							//连接节点i、节点类型j
	double P,Q,V;						//节点有功、无功（对于PV节点为无功上限）、电压（对于PV节点为需要维持电压 且为负）
};
/*****负荷*****/
struct Load_Info
{
	int i,j;							//连接节点i、节点类型j
	double P,Q,V;						//节点有功、无功、电压
};
/*****节点*****/
struct Node_Info
{
	int num;							//节点号
	int outnum;							//出线数
	int connect[20];					//与该节点相连的节点号
};
/*****PV节点*****/
struct PVNode_Info
{
	double V;							//给定电压
	int i;								//节点号
};
/*****对角元素*****/
struct Yii_Info
{
	double G,B;                         //G为实部 B为虚部
};
/*****非对角元素*****/
struct Yij_Info
{
	double G,B;                         //G为实部 B为虚部
	int j;								//非对角元素列号
};
/*****因子表上三角矩阵*****/
struct U_Info
{
	double value;
	int j;
};
/*****节点电压*****/
struct NodalVoltage_Info
{
	double V,theta;						//节点电压幅值V 角度theta
};
/*****节点功率*****/
struct NodalPower_Info
{
	double P,Q;							//节点有功P 无功Q
};


////////////////////////
////////函数声明////////
////////////////////////
void Read(char *filename);																	//数据读取
void NodeOpt();																				//节点编号优化
void YFormation();																			//导纳矩阵形成
void BFormation(Yii_Info *Yii,Yij_Info *Yij,int flag,int *NUsum,U_Info *U,double *D);		//形成B' B''及因子表
void Initial(double V0,double theta0);														//赋节点电压初值
void PQCal(int flag);																		//节点功率计算
void PQErrorCal(double *DI,int flag);																//最大功率误差计算
void SolvEqu(U_Info *U,double *D,double *DI,int *NUsum);									//解线性方程组
void OutData();






#endif