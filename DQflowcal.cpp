#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdio.h>
#include <iomanip.h>
#include <stdlib.h>
#include "definition.h"

#define Error 0.00001

////////////////////////
//////定义全局变量//////
////////////////////////

/*****数据读取*****/
int					number;						//用户输入值
char				file[10],*filename;			//用户输入文件名 文件名指针filename

/*****系统基本信息*****/
int					Mmax;						//最大迭代次数
double				SB,V0;						//基准容量SB 基准电压V0
int					N,Npv,Nb,Nt,Ng,Nl;			//节点数N PV节点数Npv 线路支路数Nb 变压器支路数Nt 发电机节点数Ng 负荷节点数Nl
int					Nbg=0;						//接地支路数
Branch_Info			*Branch;					//支路数据
Generator_Info		*Generator;					//发电机节点数据
Load_Info			*Load;						//负荷节点数据
Node_Info			*Node;						//节点数据
PVNode_Info			*PVNode;					//PV节点数据
NodalPower_Info		*NodalPower;				//节点功率
NodalPower_Info		*GenePower;					//发电机节点功率
NodalVoltage_Info	*NodalVoltage;				//节点电压

/*****节点优化*****/
int					*OptNodeNum;				//优化节点排序数组指针
int					ErrorNode;					//最大功率误差节点号

/*****导纳矩阵相关*****/
int					*NYseq,*NYsum;				//各行非对角元素首地址数组指针NYseq 各行非对角元素个数数组指针NYsum
U_Info				*U1,*U2;					//B'上三角矩阵元素U1 B''上三角矩阵元素U2
double				*D1,*D2;					//B'对角矩阵元素D1 B''对角矩阵元素D2
int					*NUsum1,*NUsum2;			//B'因子表上三角矩阵各行非对角元素数 B''因子表上三角矩阵各行非对角元素数
Yii_Info			*Yii;
Yij_Info			*Yij;						//导纳矩阵Y对角元素Yii 非对角元素Yij
Yii_Info			*Yii1;
Yij_Info			*Yij1;						//用于形成B'的矩阵Y1对角元素Yii1 非对角元素Yij1
Yii_Info			*Yii2;
Yij_Info			*Yij2;						//用于形成B''的矩阵Y2对角元素Yii2 非对角元素Yij2

double				NowError;
int					flag;
double				*DI1,*DI2;
double				Dtheta,Dv;
double				MaxError;
int					MaxErrorNode;
int					cnt;

////////////////////////
//////////函数//////////
////////////////////////

/*****主函数*****/
void main()
{
	while(1)
	{
		/*step1：显示欢迎信息*/

		cout<<endl<<"****************************PQ解耦法潮流计算程序****************************"<<endl<<endl;
		cout<<"			1.IEEE 5节点系统"<<endl;
		cout<<"			2.IEEE 14节点系统"<<endl;
		cout<<"			3.IEEE 30节点系统"<<endl;
		cout<<"			4.IEEE 57节点系统"<<endl;
		cout<<"			5.IEEE 118节点系统"<<endl;
		cout<<"			6.自定义文件"<<endl;
		cout<<"请输入序号并回车：";



		/*step2：用户输入*/

		cin>>number;
		switch(number)
		{
		case 1:
			filename="5.txt";
			break;
		case 2:
			filename="14.txt";
			break;
		case 3:
			filename="30.txt";
			break;
		case 4:
			filename="57.txt";
			break;
		case 5:
			filename="118.txt";
			break;
		case 6:
			cout<<"请输入文件名：";
			cin>>file;
			filename=file;
			break;
		}



		/*step3：确定文件名是否正确 打开文件并读取数据*/

		ifstream in(filename,ios::nocreate);
		if(!in)
		{
			cout<<"文件名错误！按任意键返回"<<endl;
			system("pause");
			continue;
		}
		else
			Read(filename);



		/*step3：T2节点编号优化*/

		NodeOpt();



		/*step4：形成导纳矩阵*/

		YFormation();



		/*step5：形成B' B''及因子表*/

		NUsum1=new int[N+1];						
		NUsum2=new int[N+1];
		D1=new double[N+1];
		D2=new double[N+1];
		U1=new U_Info[(Nb+Nt)*8];
		U2=new U_Info[(Nb+Nt)*8];
		BFormation(Yii1,Yij1,1,NUsum1,U1,D1);
		BFormation(Yii2,Yij2,2,NUsum2,U2,D2);

		//输出
		cout<<endl<<"5.生成因子表"<<endl;
		//输出上三角矩阵各行非对角元素数
		cout<<endl<<"第一因子表上三角矩阵各行非对角元素数："<<endl;
		for(int i=1;i<N;i++)
			cout<<"第"<<i<<"行	"<<NUsum1[i]<<"个"<<endl;
		cout<<endl<<"第二因子表上三角矩阵各行非对角元素数："<<endl;
		for(i=1;i<N;i++)
			cout<<"第"<<i<<"行	"<<NUsum2[i]<<"个"<<endl;
		cout<<endl;
		system("pause");
		//输出因子表对角元素
		cout<<endl<<"因子表对角元素D1:"<<endl;
		int t1=0;
		int t2=0;
		for(i=1;i<N;i++)
		{
			cout<<i<<"	"<<D1[i]<<endl;
			if(i<N)
			{	
				t1+=NUsum1[i];
				t2+=NUsum2[i];
			}
		}
		cout<<endl<<"因子表对角元素D2:"<<endl;
		for(i=1;i<N;i++)
			cout<<i<<"	"<<D2[i]<<endl;
		cout<<endl;
		system("pause");
		//输出因子表上三角元
		cout<<endl<<"第一因子表上三角元U1="<<endl;
		for(i=1;i<=t1;i++)
			cout<<i<<setw(12)<<U1[i].value<<endl;
		cout<<endl<<"第二因子表上三角元U2="<<endl;
		for(i=1;i<=t2;i++)
			cout<<i<<setw(12)<<U2[i].value<<endl;
		cout<<endl;
		system("pause");



		/*step6：赋节点电压初值*/

		NodalVoltage=new NodalVoltage_Info[N+1];
		Initial(1,0);



		/*step7：进行迭代*/
		cnt=1;
		NodalPower=new NodalPower_Info[N+1];
		GenePower=new NodalPower_Info[N+1];
		DI1=new double[N+1];
		DI2=new double[N+1];
		NowError=100;									//赋误差初值
		cout<<endl<<"按任意键开始迭代"<<endl<<endl;
		system("pause");
		while(NowError>Error)							//误差不满足要求则继续迭代
		{
			//P迭代
			PQCal(1);									//计算节点功率
			PQErrorCal(DI1,1);							//计算节点功率误差
			SolvEqu(U1,D1,DI1,NUsum1);					//解线性方程组
			for(i=1;i<N;i++)							//修正
			{
				Dtheta=DI1[i]/NodalVoltage[i].V;
				NodalVoltage[i].theta-=Dtheta;
			}
			NowError=MaxError;
			MaxErrorNode=ErrorNode;

			//Q迭代
			PQCal(2);									//计算节点功率
			PQErrorCal(DI2,2);							//计算节点功率误差
			SolvEqu(U2,D2,DI2,NUsum2);					//解线性方程组
			for(i=1;i<N;i++)							//修正
			{
				Dv=DI2[i];
				NodalVoltage[i].V-=Dv;
			}

			if(NowError<MaxError)						//判断最大误差
			{
				NowError=MaxError;
				MaxErrorNode=ErrorNode;
			}

			//输出每次迭代结果
			cout<<"============================================================="<<endl;
			cout<<endl<<"第"<<cnt<<"次迭代"<<endl;
			cout<<"最大误差："<<NowError<<"	最大误差节点号："<<MaxErrorNode<<endl<<endl;
			for(i=1;i<N;i++)
				cout<<"节点"<<i<<"	ΔP:"<<setw(15)<<DI1[i]<<setw(5)<<"θ:"<<setw(15)<<NodalVoltage[i].theta<<endl;
			for(i=1;i<N;i++)
				cout<<"节点"<<i<<"	ΔQ:"<<setw(15)<<DI2[i]<<setw(5)<<"V:"<<setw(15)<<NodalVoltage[i].V<<endl;
			cnt++;
			cout<<endl;
			if(NowError<Error)
			{
				cout<<"经"<<cnt<<"次迭代后收敛"<<endl<<endl;
			}
			if(cnt>Mmax)
			{
				cout<<"潮流计算结果不收敛!"<<endl;
				break;
			}
		}



		/*step8：输出结果*/
		OutData();
		system("pause");
	}
}



/*****数据读取*****/
void Read(char *filename)
{
	char head[100],enter;						//表头 换行符寄存器



	/*step1：打开文件*/

	ifstream in(filename);



	/*step2：读取信息*/

	//基本信息
	in.getline(head,100);						//跳过表头
	in>>Mmax>>SB>>V0;							//读取最大迭代次数Mmax 系统基准容量SB 系统基准电压V0
	in>>enter;									//跳过换行符
	in.getline(head,100);						//跳过表头
	in>>N>>Nb>>Nt>>Ng>>Nl;						//读取节点数N 线路支路数Nb 变压器支路数Nt 发电机节点数Ng 负荷节点数Nl

	//支路
	in>>enter;									//跳过换行符
	in.getline(head,100);						//跳过表头
	Branch=new Branch_Info[Nb+Nt+1];
	for(int i=1;i<=Nb+Nt;i++)
	{
		in>>Branch[i].i>>Branch[i].j>>Branch[i].R>>Branch[i].X>>Branch[i].YK;
	}

	//发电机节点
	in>>enter;									//跳过换行符
	in.getline(head,100);						//跳过表头
	Generator=new Generator_Info[Ng+1];
	for(i=1;i<=Ng;i++)
	{
		in>>Generator[i].i>>Generator[i].j>>Generator[i].P>>Generator[i].Q>>Generator[i].V;
	}

	//负荷节点
	Load=new Load_Info[Nl+1];
	for(i=1;i<=Nl;i++)
	{
		in>>Load[i].i>>Load[i].j>>Load[i].P>>Load[i].Q>>Load[i].V;
	}
    in.close();



	/*step3:输出*/

	cout<<endl<<"***************************************************************"<<endl;
	cout<<endl<<"1.系统信息"<<endl;

	//基本信息
	cout<<endl<<"=======================基本信息======================="<<endl;
	cout<<"最大迭代次数\t"<<"系统基准容量\t"<<"系统基准电压"<<endl;
	cout<<Mmax<<"\t\t"<<SB<<"\t\t"<<V0<<endl;
	cout<<"系统节点数\t"<<"线路支路数\t"<<"变压器支路数\t"<<"发电机节点数\t"<<"负荷节点数"<<endl;
	cout<<N<<"\t\t"<<Nb<<"\t\t"<<Nt<<"\t\t"<<Ng<<"\t\t"<<Nl<<endl<<endl;
	system("pause");

	//支路信息
	cout<<endl<<"=======================支路信息======================="<<endl;
	cout<<"支路连接点i\t"<<"支路连接点j\t"<<"支路电阻\t"<<"支路电抗\t"<<"支路接地导纳"<<endl;
	for(i=1;i<=Nb+Nt;i++)
	{
		cout<<Branch[i].i<<"\t\t"<<Branch[i].j<<"\t\t"<<Branch[i].R<<"\t\t"<<Branch[i].X<<"\t\t"<<Branch[i].YK<<endl;
	}
	cout<<endl;
	system("pause");

	/*节点信息*/
	cout<<endl<<"=======================节点信息======================="<<endl;
	cout<<"节点编号\t"<<"节点类型\t"<<"节点有功功率\t"<<"节点无功功率\t"<<"节点电压"<<endl;
	for(i=1;i<=Ng;i++)
	{
		cout<<Generator[i].i<<"\t\t"<<Generator[i].j<<"\t\t"<<Generator[i].P<<"\t\t"<<Generator[i].Q<<"\t\t"<<Generator[i].V<<endl;
	}
	for(i=1;i<=Nl;i++)
	{
		cout<<Load[i].i<<"\t\t"<<Load[i].j<<"\t\t"<<Load[i].P<<"\t\t"<<Load[i].Q<<"\t\t"<<Load[i].V<<endl;
	}
	cout<<endl;
	system("pause");
}

/*****节点编号优化*****/
void NodeOpt()
{
	Node=new Node_Info[N+1];					//节点结构体
	int Balance=0;								//平衡节点号



	/*step1：初始化*/

	for(int i=1;i<=N;i++)
	{
		Node[i].num=i;
		Node[i].outnum=0;
		for(int ii=1;ii<=20;ii++)
			Node[i].connect[ii]=0;
	}



	/*step2：确定平衡节点 且将其出线数设为1000*/

	for(i=1;i<=Ng;i++)
	{
		if(Generator[i].j==0)
			Balance=Generator[i].i;
	}
	Node[Balance].outnum=1000;



	/*step3：确定各节点连接关系*/

	for(i=1;i<=Nb+Nt;i++)            
	{
		//判断是否为接地支路 若是则不添加连接关系
		if(Branch[i].i==Branch[i].j)
		{
			Nbg++;
			continue;
		}

		//判断是否为连接平衡节点的支路 若是则不添加连接关系
		else if(Branch[i].i==Balance||Branch[i].i==Balance)		
			continue;

		//判断是否已有连接关系 若是则不再添加连接关系
		for(int ii=1;ii<=Node[Branch[i].i].outnum;ii++)
		{
			if(Node[Branch[i].i].connect[ii]==Branch[i].j)
				break;
		}

		//添加连接关系
		if(ii>Node[Branch[i].i].outnum)							
		{
			Node[Branch[i].i].outnum++;
			Node[Branch[i].i].connect[Node[Branch[i].i].outnum]=Branch[i].j;
			Node[Branch[i].j].outnum++;
			Node[Branch[i].j].connect[Node[Branch[i].j].outnum]=Branch[i].i;
		}
	}



	/*step4：输入PV节点信息*/

	Npv=0;
	PVNode=new PVNode_Info[Ng+1];
	for(i=1;i<=Ng;i++)
	{
		if(Generator[i].j==-1)
		{
			Npv++;
			PVNode[Npv].i=Generator[i].i;
			PVNode[Npv].V=Generator[i].V;
		}
	}



	/*step5：输出各节点连接情况*/

	cout<<endl<<"***************************************************************"<<endl;
	cout<<endl<<"2.节点连接情况"<<endl<<endl;

	for(i=1;i<=N;i++)
	{
		if(Node[i].outnum>0&&i!=Balance)
		{
			cout<<"节点"<<i<<"出线数为"<<Node[i].outnum<<"，连接的节点有：";
			for(int ii=1;ii<=Node[i].outnum;ii++)
				cout<<Node[i].connect[ii]<<" ";
			cout<<endl;
		}
	}
	cout<<endl<<"PV节点有"<<Npv<<"个，为";
	for(i=1;i<=Npv;i++)
	{
		cout<<PVNode[i].i<<" ";
	}
	cout<<endl<<"平衡节点为"<<Balance<<endl;
	cout<<"接地支路有"<<Nbg<<"条"<<endl<<endl;
	system("pause");



	/*step6：半动态节点编号优化*/

	OptNodeNum=new int[N+1];
	int leastnum;
	OptNodeNum[N]=Balance;
	for(i=1;i<=N-1;i++)
	{
		//找出出线数最小节点
		leastnum=1;
		for(int j=1;j<=N;j++)								
		{
			if(Node[j].outnum<Node[leastnum].outnum)
				leastnum=j;
		}
		OptNodeNum[i]=leastnum;												//将新节点序号写入优化编号数组
		/*
		cout<<"此次消去节点为："<<leastnum<<endl;							//输出消去节点
		system("pause");
		*/

		//消去节点 更新出线数
		for(j=1;j<=Node[leastnum].outnum;j++)
		{
			for(int k=1;k<=Node[Node[leastnum].connect[j]].outnum;k++)
				if(Node[Node[leastnum].connect[j]].connect[k]==leastnum)	//判断消去节点是否是该节点连接列表的最后一个
					break;
			if(k==Node[Node[leastnum].connect[j]].outnum)					//若消去节点是该节点连接列表的最后一个 则接将该节点出线数-1
				Node[Node[leastnum].connect[j]].outnum--;
			else															//若消去节点不是该节点连接列表的最后一个 则将该节点出线数-1 并将其他相连节点前移
			{
			for(k=k;k<Node[Node[leastnum].connect[j]].outnum;k++)
				Node[Node[leastnum].connect[j]].connect[k]=Node[Node[leastnum].connect[j]].connect[k+1];
			Node[Node[leastnum].connect[j]].outnum--;
			}
		}
		/*
		for(j=1;j<=N;j++)												//输出消去后节点连接变化情况
		{
			if(Node[j].outnum!=1000&&j!=leastnum)
			{
				cout<<"此次消去后节点"<<j<<"的出线数为"<<Node[j].outnum<<"，与其连接的节点有：";
				for(int k=1;k<=Node[j].outnum;k++)
				cout<<"节点"<<Node[j].connect[k]<<" ";
			cout<<endl;
			}
		}
		system("pause");
		*/

		//追加支路 更新出线数
		for(j=1;j<=Node[leastnum].outnum-1;j++)							//判断与消去节点leastnum相连的第j节点和第k节点是否相连 若是则追加支路
		{
			for(int k=j+1;k<=Node[leastnum].outnum;k++)
			{
				for(int n=1;n<=Node[Node[leastnum].connect[j]].outnum;n++)
				{
					if(Node[Node[leastnum].connect[j]].connect[n]==Node[leastnum].connect[k])
						break;
				}
				if(n>Node[Node[leastnum].connect[j]].outnum)
				{
					Node[Node[leastnum].connect[j]].outnum++;
					Node[Node[leastnum].connect[j]].connect[Node[Node[leastnum].connect[j]].outnum]=Node[leastnum].connect[k];
					Node[Node[leastnum].connect[k]].outnum++;
					Node[Node[leastnum].connect[k]].connect[Node[Node[leastnum].connect[k]].outnum]=Node[leastnum].connect[j];
				}
			}
		}
		/*
		for(j=1;j<=N;j++)												//输出追加后节点连接变化情况
		{
			if(Node[j].outnum!=1000&&j!=leastnum)
			{
				cout<<"此次追加后节点"<<j<<"的出线数为"<<Node[j].outnum<<"，与其连接的节点有：";
				for(int k=1;k<=Node[j].outnum;k++)
					cout<<"节点"<<Node[j].connect[k]<<" ";
				cout<<endl;
			}
		}
		system("pause");
		*/
		Node[leastnum].outnum=1000;										//将已排序过节点排除比较
	}



	/*step7：更新节点编号*/

	for(i=1;i<=N;i++)
		Node[OptNodeNum[i]].num=i;

	//PV节点
	for(i=1;i<=Npv;i++)
		PVNode[i].i=Node[PVNode[i].i].num;
	//发电机节点
	for(i=1;i<=Ng;i++)
		Generator[i].i=Node[Generator[i].i].num;
	//负荷节点
	for(i=1;i<=Nl;i++)
		Load[i].i=Node[Load[i].i].num;


	/*step8：节点重新排序*/
	//PV节点
	PVNode_Info PVtemp;
	for(i=1;i<Npv;i++)
		for(int j=Npv;j>i;j--)
			if(PVNode[j-1].i>PVNode[j].i)
			{
				PVtemp=PVNode[j];
				PVNode[j]=PVNode[j-1];
				PVNode[j-1]=PVtemp;
			}
	//发电机节点
	Generator_Info Genetemp;
	for(i=1;i<Ng;i++)
		for(int j=Ng;j>i;j--)
			if(Generator[j-1].i>Generator[j].i)
			{
				Genetemp=Generator[j];
				Generator[j]=Generator[j-1];
				Generator[j-1]=Genetemp;
			}
	//负荷节点
	Load_Info Loadtemp;
	for(i=1;i<Nl;i++)
		for(int j=Nl;j>i;j--)
			if(Load[j-1].i>Load[j].i)
			{
				Loadtemp=Load[j];
				Load[j]=Load[j-1];
				Load[j-1]=Loadtemp;
			}



	/*step9：更新支路和节点信息*/

	//线路支路
	for(i=1;i<=Nb;i++)													
	{
		Branch[i].i=Node[Branch[i].i].num;
		Branch[i].j=Node[Branch[i].j].num;
	}
	//变压器支路
	for(i=Nb+1;i<=Nb+Nt;i++)
	{
		Branch[i].i=Node[Branch[i].i].num;
		Branch[i].j=-Node[Branch[i].j].num;								//负号作为变压器支路标志
	}
	//支路编号更新
	Branch_Info Branchtem;												//临时变量
	int k=0;
	for(i=1;i<=Nb+Nt-Nbg;i++)         
	{
		if(abs(Branch[i].i)>abs(Branch[i].j))							//支路两端节点号较小的排在前面
		{
			int j=Branch[i].i;
			Branch[i].i=Branch[i].j;
			Branch[i].j=j;

		}
		else if(Branch[i].i==Branch[i].j)								//将接地支路排在最后
		{
			Branchtem=Branch[i];
			Branch[i]=Branch[Nb+Nt-k];
			Branch[Nb+Nt-k]=Branchtem;
			k++;
			if(abs(Branch[i].i)>abs(Branch[i].j))   
			{
				int j=Branch[i].i;
				Branch[i].i=Branch[i].j;
				Branch[i].j=j;
			}
		}
	}
	k=1;
	for(i=1;i<=N;i++)													//将支路按照节点编号从小到大排列
	{
		for(int j=1;j<=Nb+Nt-Nbg;j++)
		{
			if(abs(Branch[j].i)==i)
			{
				Branchtem=Branch[k];
				Branch[k]=Branch[j];
				Branch[j]=Branchtem;
				k++;
			}
		}
	}
	for(i=1;i<=N;i++)
	{
		for(int j=i+1;j<=Nb+Nt-Nbg;j++)
		{
			if(abs(Branch[i].i)!=abs(Branch[j].i))
				break;
			if(abs(Branch[i].i)==abs(Branch[j].i)&&abs(Branch[i].j)>abs(Branch[j].j))
			{
				Branchtem=Branch[i];
				Branch[i]=Branch[j];
				Branch[j]=Branchtem;
			}
		}
	}



	/*step10：输出优化结果*/

	cout<<endl<<"***************************************************************"<<endl;
	cout<<endl<<"3.T2节点编号优化"<<endl<<endl;

	//节点优化结果
	for(i=1;i<=N;i++)
		cout<<"原第"<<i<<"节点，经半动态优化后编号为"<<Node[i].num<<endl;
	cout<<endl;
	system("pause");
	//优化后支路信息
	cout<<endl<<"=======================优化后支路信息======================="<<endl;
	cout<<"支路连接点i\t"<<"支路连接点j\t"<<"支路电阻\t"<<"支路电抗\t"<<"支路接地导纳"<<endl;
	for(i=1;i<=Nb+Nt;i++)
	{
		cout<<Branch[i].i<<"\t\t"<<Branch[i].j<<"\t\t"<<Branch[i].R<<"\t\t"<<Branch[i].X<<"\t\t"<<Branch[i].YK<<endl;
	}
	cout<<endl;
	system("pause");
	//优化后节点信息
	cout<<endl<<"=======================优化后节点信息======================="<<endl;
	cout<<"节点编号\t"<<"节点类型\t"<<"节点有功功率\t"<<"节点无功功率\t"<<"节点电压"<<endl;
	for(i=1;i<=Ng;i++)
	{
		cout<<Generator[i].i<<"\t\t"<<Generator[i].j<<"\t\t"<<Generator[i].P<<"\t\t"<<Generator[i].Q<<"\t\t"<<Generator[i].V<<endl;
	}
	for(i=1;i<=Nl;i++)
	{
		cout<<Load[i].i<<"\t\t"<<Load[i].j<<"\t\t"<<Load[i].P<<"\t\t"<<Load[i].Q<<"\t\t"<<Load[i].V<<endl;
	}
	cout<<endl;
	cout<<"PV节点有"<<Npv<<"个，为";
	for(i=1;i<=Npv;i++)
	{
		cout<<PVNode[i].i<<" ";
	}
	cout<<endl;
	cout<<"平衡节点为节点"<<N<<endl;
	cout<<"接地支路有"<<Nbg<<"条"<<endl<<endl;
	system("pause");
}

/*****导纳矩阵形成*****/
void YFormation()
{
	int i,j;															//连接节点编号
	double R,X,YK,Zmag2;												//支路电阻 电抗 接地导纳 阻抗幅值
	double G,B,b;														//导纳实部 虚部
	Yii=new Yii_Info[N+1];
	Yij=new Yij_Info[Nb+Nt+1];											//Y存放导纳矩阵
	Yii1=new Yii_Info[N+1];
	Yij1=new Yij_Info[Nb+Nt+1];											//Y1用于形成B'（不考虑接地支路）
	Yii2=new Yii_Info[N+1];
	Yij2=new Yij_Info[Nb+Nt+1];											//Y2用于形成B''（不考虑线路电阻）
	NYsum=new int[N+1];
	NYseq=new int[N+1];



	/*step1：初始化*/
	for(int n=1;n<=N;n++)												//初始化对角元素
	{
		Yii[n].G=0;
		Yii[n].B=0;
		Yii1[n].G=0;
		Yii1[n].B=0;
		Yii2[n].G=0;
		Yii2[n].B=0;
		NYsum[n]=0;
	}
	for(n=1;n<=Nb+Nt;n++)											//初始化非对角元素
	{
		Yij[n].G=0;
		Yij[n].B=0;
		Yij1[n].G=0;
		Yij1[n].B=0;
		Yij2[n].G=0;
		Yij2[n].B=0;
	}



	/*step2：形成不考虑接地支路的导纳矩阵*/
	for(n=1;n<=Nb+Nt;n++)
	{
		i=abs(Branch[n].i);
		j=abs(Branch[n].j);
		R=Branch[n].R;
		X=Branch[n].X;
		YK=Branch[n].YK;

		if(i==j)														//跳过对地支路
			continue;

		//计算该支路导纳
		Zmag2=R*R+X*X;													//阻抗幅值
		G=R/Zmag2;														//导纳实部
		B=-X/Zmag2;														//导纳虚部
		b=-1/X;															//不考虑电阻的导纳虚部

		if(Branch[n].i<0||Branch[n].j<0)
		{
			Yij[n].G=-G/YK;												//计算变压器支路互导纳 即非对角元素
			Yij[n].B=-B/YK;
			Yij1[n].G=-G/YK;
			Yij1[n].B=-B/YK;
			Yij2[n].G=0;
			Yij2[n].B=-b/YK;

			Yii[i].G=Yii[i].G+G/YK;										//计算变压器节点自导纳 即对角元素
			Yii[i].B=Yii[i].B+B/YK;
			Yii[j].G=Yii[j].G+G/YK;
			Yii[j].B=Yii[j].B+B/YK;
			Yii1[i].G=Yii1[i].G+G/YK;
			Yii1[i].B=Yii1[i].B+B/YK;
			Yii1[j].G=Yii1[j].G+G/YK;
			Yii1[j].B=Yii1[j].B+B/YK;
			Yii2[i].B=Yii2[i].B+b/YK;
			Yii2[j].B=Yii2[j].B+b/YK;
		}
		else
		{
			Yij[n].G=-G;												//计算线路支路互导纳 即非对角元素
			Yij[n].B=-B;
			Yij1[n].G=-G;
			Yij1[n].B=-B;
			Yij2[n].G=0;
			Yij2[n].B=-b;

			Yii[i].G=Yii[i].G+G;										//计算线路节点自导纳 即对角元素
			Yii[i].B=Yii[i].B+B;
			Yii[j].G=Yii[j].G+G;
			Yii[j].B=Yii[j].B+B;
			Yii1[i].G=Yii1[i].G+G;
			Yii1[i].B=Yii1[i].B+B;
			Yii1[j].G=Yii1[j].G+G;
			Yii1[j].B=Yii1[j].B+B;
			Yii2[i].B=Yii2[i].B+b;
			Yii2[j].B=Yii2[j].B+b;
		}


		Yij[n].j=j;
		Yij1[n].j=j;
		Yij2[n].j=j;
		NYsum[i]=NYsum[i]+1;
	}
	NYseq[1]=1;
	for(i=1;i<N;i++)
		NYseq[i+1]=NYseq[i]+NYsum[i];



	/*step3：追加接地支路*/
	for(n=1;n<=Nb+Nt;n++)
	{
		i=Branch[n].i;
		j=Branch[n].j;
		R=Branch[n].R;
		X=Branch[n].X;
		YK=Branch[n].YK;

		if(i==j)													//加入对地支路
		{
			Yii[i].B=Yii[i].B+1.0/X;
			Yii2[i].B=Yii2[i].B+1.0/X;
			continue;
		}

		if(i<0||j<0)
		{
			if(i<0)
			{
				i=abs(i);											//i为非标准变比侧
				G=Yij[n].G;
				B=Yij[n].B;
				b=Yij2[n].B;

				Yii[i].G=Yii[i].G+(1.0-1.0/YK)*G;
				Yii[i].B=Yii[i].B+(1.0-1.0/YK)*B;
				Yii2[i].B=Yii2[i].B+(1.0-1.0/YK)*b;

				Yii[j].G=Yii[j].G+(1.0-YK)*G;
				Yii[j].B=Yii[j].B+(1.0-YK)*B;
				Yii2[j].B=Yii2[j].B+(1.0-YK)*b;
			}
			else
			{
				j=abs(j);											//j为非标准变比侧
				G=Yij[n].G;
				B=Yij[n].B;
				b=Yij2[n].B;

				Yii[j].G=Yii[j].G+(1.0-1.0/YK)*G;
				Yii[j].B=Yii[j].B+(1.0-1.0/YK)*B;
				Yii2[j].B=Yii2[j].B+(1.0-1.0/YK)*b;

				Yii[i].G=Yii[i].G+(1.0-YK)*G;
				Yii[i].B=Yii[i].B+(1.0-YK)*B;
				Yii2[i].B=Yii2[i].B+(1.0-YK)*b;
			}
		}
		else
		{
			B=YK;													//线路
			b=YK;
			Yii[i].B=Yii[i].B+B;
			Yii[j].B=Yii[j].B+B;
			Yii2[i].B=Yii2[i].B+b;
			Yii2[j].B=Yii2[j].B+b;
		}
	}



	/*step4：输出导纳矩阵*/
	cout<<endl<<"4.计算导纳矩阵"<<endl;
	//输出Y
	cout<<endl<<"稀疏导纳矩阵对角元:Yii="<<endl;
	for(i=1;i<=N;i++)
		cout<<i<<setw(10)<<Yii[i].G<<"+"<<setw(10)<<Yii[i].B<<"j"<<endl;
	cout<<"稀疏导纳矩阵非对角元:Yij="<<endl;
	for(i=1;i<=Nb+Nt-Nbg;i++)
	{
		cout<<i<<setw(10)<<Yij[i].G<<"+"<<setw(10)<<Yij[i].B<<"j"<<endl;
	}
	cout<<endl;
	system("pause");
	//输出Y1
	cout<<endl<<"稀疏导纳矩阵对角元:Yii1="<<endl;
	for(i=1;i<=N;i++)
		cout<<i<<setw(10)<<Yii1[i].G<<"+"<<setw(10)<<Yii1[i].B<<"j"<<endl;
	cout<<"稀疏导纳矩阵非对角元:Yij1="<<endl;
	for(i=1;i<=Nb+Nt-Nbg;i++)
	{
		cout<<i<<setw(10)<<Yij1[i].G<<"+"<<setw(10)<<Yij1[i].B<<"j"<<endl;
	}
	cout<<endl;
	system("pause");
	//输出Y2
	cout<<endl<<"稀疏导纳矩阵对角元:Yii2="<<endl;
	for(i=1;i<=N;i++)
		cout<<i<<setw(10)<<Yii2[i].G<<"+"<<setw(10)<<Yii2[i].B<<"j"<<endl;
	cout<<"稀疏导纳矩阵非对角元:Yij="<<endl;
	for(i=1;i<=Nb+Nt-Nbg;i++)
	{
		cout<<i<<setw(10)<<Yij2[i].G<<"+"<<setw(10)<<Yij2[i].B<<"j"<<endl;
	}
	cout<<endl;
	system("pause");
}

/*****形成B' B''及因子表*****/
void BFormation(Yii_Info *Yii,Yij_Info *Yij,int flag,int *NUsum,U_Info *U,double *D)
//矩阵对角元素Yii 非对角元素Yij Yii1和Yij1用于形成B' Yii2和Yij2用于形成B''
//标识变量flag=1形成B'因子表 flag=2形成B''因子表
//NUsum存放上三角矩阵各行非对角元素数
//D存放因子表的对角元素
{
	int n_pv,i_pv;													//节点数组计数n_pv 节点号i_pv
	int i,j;														//行号变量i 列下标变量j
	int n,n_u;														//临时计数变量n 因子表上三角矩阵元素计数变量
	int i_above;													//消去行号i计数变量
	double	*B;														//系数矩阵B元素
	double Btemp;
	B=new double[N+1];												//临时变量
	n_pv=1;
	i_pv=PVNode[1].i;

	for(i=1;i<N;i++)
	{
		if(flag==2&&i==i_pv)
		{
			//F部分
			n_pv++;
			i_pv=PVNode[n_pv].i;
			NUsum[i]=0;
			D[i]=0.0;
		}
		else
		{
			//A部分
			for(n=i+1;n<N;n++)
				B[n]=0;
			B[i]=Yii[i].B;

			//B部分
			for(n=NYseq[i];n<NYseq[i+1];n++)
			{
				j=Yij[n].j;
				B[j]=Yij[n].B;
			}

			//C部分
			if(flag==2)
			{
				for(n=1;n<=Npv;n++)
					B[PVNode[n].i]=0;
			}

			//D部分
			n_u=1;
			for(i_above=1;i_above<i;i_above++)
			{
				n=1;
				while(n<=NUsum[i_above])
				{
					if(U[n_u].j==i)
					{
						Btemp=U[n_u].value/D[i_above];
						while(n<=NUsum[i_above])
						{
							j=U[n_u].j;
							B[j]-=Btemp*U[n_u].value;
							n++;
							n_u++;
						}
					}
					else
					{
						n++;
						n_u++;
					}
				}
			}

			//E部分
			Btemp=1.0/B[i];
			D[i]=Btemp;
			n=0;
			for(j=i+1;j<N;j++)
			{
				if(B[j]!=0)
				{
					U[n_u].value=B[j]*Btemp;
					U[n_u].j=j;
					n++;
					n_u++;
				}
			}
			NUsum[i]=n;
		}
	}
}



/*赋节点电压初值*/
void Initial(double V0,double theta0)
{
	for(int i=1;i<=N;i++)
	{
		NodalVoltage[i].V=V0;
		NodalVoltage[i].theta=theta0;
	}
	for(i=1;i<=Npv;i++)
	{   
		NodalVoltage[PVNode[i].i].V=PVNode[i].V;
	}
}



/*节点功率计算*/
void PQCal(int flag)							//flag=1进行P-theta迭代 flag=2进行Q-V迭代
{
	double Vi,A,B,vv,theta;						//i节点电压幅值Vi 临时变量A B vv theta
	int j;										//临时变量j
	for(int i=1;i<=N;i++)
	{
		if (flag==1)
			NodalPower[i].P=0;
		else
			NodalPower[i].Q=0;
	}
	for(i=1;i<=N;i++)
	{
		Vi=NodalVoltage[i].V;
		if(flag==1)
		{
			A=Yii[i].G;                   
			NodalPower[i].P=NodalPower[i].P+Vi*Vi*A;
		}
		else
		{
			A=-Yii[i].B;
			NodalPower[i].Q=NodalPower[i].Q+Vi*Vi*A;
		}
		if(i==N)								//最后一行无非对角元，跳出
			break;
		for(int n=NYseq[i];n<NYseq[i+1];n++)	//找到第i行的非对角元素
		{
			A=Yij[n].G;
			B=Yij[n].B;
			j=Yij[n].j;
			vv=Vi*NodalVoltage[j].V;
			theta=NodalVoltage[i].theta-NodalVoltage[j].theta;
			if(flag==1)
			{
				A=A*vv*cos(theta);
				B=B*vv*sin(theta);
				NodalPower[i].P=NodalPower[i].P+A+B;
				NodalPower[j].P=NodalPower[j].P+A-B;
			}
			else
			{
				B=B*vv*cos(theta);
				A=A*vv*sin(theta);
				NodalPower[i].Q=NodalPower[i].Q+A-B;
				NodalPower[j].Q=NodalPower[j].Q-A-B;
			}
		}
	}
}

/*最大功率误差计算*/
void PQErrorCal(double *DI,int flag)								//功率误差DI flag=1计算P误差 flag=2计算Q误差
{
	MaxError=0;												//最大功率误差MaxError
	double Vi;														//i点电压幅值Vi
	double Wi,Wtemp;
	int i=1;
	int n_g=1;
	int n_l=1;
	int n_pv=1;
	int i_g=Generator[1].i;
	int i_l=Load[1].i;
	int i_pv=PVNode[1].i;
	for(i=1;i<=N;i++)
	{
		Vi=NodalVoltage[i].V;

		//负荷
		if(i==i_l)
		{
			if(flag==1)
				Wi=-Load[n_l].P;
			else
				Wi=-Load[n_l].Q;
			n_l++;
			i_l=Load[n_l].i;
		}
		else
			Wi=0;
		Wtemp=Wi;

		//节点注入功率
		if(flag==1)
			Wi=Wi-NodalPower[i].P;
		else
			Wi=Wi-NodalPower[i].Q;

		//发电机
		if(i==i_g)
		{
			if(flag==1)
			{
				NodalPower[i].P=Wtemp;
				GenePower[i].P=-Wi;
				Wi=Wi+Generator[n_g].P;
			}
			else
			{
				NodalPower[i].Q=Wtemp;
				GenePower[i].Q=-Wi;
				Wi=Wi+Generator[n_g].Q;
			}
			n_g++;
			i_g=Generator[n_g].i;
		}

		if(i==N)									//平衡节点不计算功率误差
			break;

		if(flag==2&&i==i_pv)
		{
			n_pv++;
			i_pv=PVNode[n_pv].i;
			DI[i]=0;
		}
		else
		{
			if(fabs(Wi)>MaxError)
			{
				MaxError=fabs(Wi);
				ErrorNode=i;
			}
			DI[i]=Wi/Vi;
		}
	}
}

/*解线性方程组*/
void SolvEqu(U_Info *U,double *D,double *DI,int *NUsum)
{
	int	n_u=1;
	int i,j,count;
	double	DItemp;
	for(i=1;i<N;i++)
	{
		DItemp=DI[i];
		for(count=1;count<=NUsum[i];count++)				//前代
		{
			j=U[n_u].j;
			DI[j]=DI[j]-DItemp*U[n_u].value;
			n_u++;
		}
		DI[i]=DItemp*D[i];
	}
	for(i=N-1;i>0;i--)										//回代
	{
		DItemp=DI[i];
		for(count=1;count<=NUsum[i];count++)         
		{
			n_u--;
			j=U[n_u].j;
			DItemp=DItemp-DI[j]*U[n_u].value;
		}
		DI[i]=DItemp;
	}
}

///////////////////////////数据输出子函数////////////////////////////////////////
void OutData()
{
	ofstream	output;											//输出
	int		i,j,k,n,i_g,n_g,VminNode;
	double	theta,Vmin,V,P,Q;
	double	PLoss,QLoss,R,X,YK,Vi,Vj,Ei,Ej,Fi,Fj,DE,DF,Ir,Ii,Pij,Qij,Pji,Qji,Zmag;
	Vmin=NodalVoltage[1].V;           //系统最低电压      
	i_g=Generator[1].i;       
	VminNode=1;					//系统最低电压 所在节点号                     
	n_g=1;
	output.open("output.txt");
	///////////////////////////输出///////////////////////////////
	output<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(6); //屏幕输出格式控制。
	output<<"******************************节点信息输出：*******************************"<<endl;
	output<<setw(4)<<"节点号"<<setw(12)<<"节点电压V"<<setw(12)<<"相角θ"<<setw(12)<<"发电机有功"<<setw(12)<<"发电机无功"<<setw(12)<<"PL"<<setw(12)<<"QL"<<endl;
	for(i=1;i<N+1;i++)
	{
		k=OptNodeNum[i];
		theta=NodalVoltage[i].theta*180/3.14;
		V=NodalVoltage[i].V;
		if(V<Vmin)
		{
			Vmin=V;
			VminNode=i;
		}
 		if(i==i_g)
		{
			P=GenePower[i].P;
			Q=GenePower[i].Q;
			n_g=n_g+1;
			i_g=Generator[n_g].i;
		}
		else
		{
			P=0;
			Q=0;
		}
		output<<setw(4)<<k<<setw(12)<<V<<setw(12)<<theta<<setw(12)<<P<<setw(12)<<Q<<setw(12)<<NodalPower[k].P<<setw(12)<<NodalPower[k].Q<<endl;
	}
	PLoss=0;
	QLoss=0;
	output<<endl<<"**********************支路数据输出：**************************"<<endl;
	output<<"   i"<<"   j"<<"     Pij"<<"         Qij"<<"          Pji"<<"         Qji"<<endl;
	for(n=1;n<=Nb+Nt;n++)
	{
		i=abs(Branch[n].i);
		j=abs(Branch[n].j);
		R=Branch[n].R;
		X=Branch[n].X;
		YK=Branch[n].YK;
		Vi=NodalVoltage[i].V;
		theta=NodalVoltage[i].theta;
		Ei=Vi*cos(theta);
		Fi=Vi*sin(theta);
		Vj=NodalVoltage[j].V;
		theta=NodalVoltage[j].theta;
		Ej=Vj*cos(theta);
		Fj=Vj*sin(theta);
		if(Branch[n].i<0||Branch[n].j<0)
		{
			if(Branch[n].i<0)
			{
				Ei=Ei/YK;
				Fi=Fi/YK;
			}
			else
			{
				Ej=Ej/YK;
				Fi=Fj/YK;
			}
			YK=0;
		}
		DE=Ei-Ej;
		DF=Fi-Fj;
		Zmag=R*R+X*X;
		Ir=(DE*R+DF*X)/Zmag;
		Ii=(DF*R-DE*X)/Zmag;
		Pij=Ir*Ei+Ii*Fi;
		Qij=Ir*Fi-Ii*Ei;
		Pji=-Ir*Ej-Ii*Fj;
		Qji=-Ir*Fj+Ii*Ej;
		Qij-=Vi*Vi*YK/2;
		Qji-=Vj*Vj*YK/2;
		PLoss+=Pij+Pji;
		QLoss+=Qij+Qji;
		i=OptNodeNum[i];
		j=OptNodeNum[j];
		if(i>j)
		{
			k=j;
			j=i;
			i=k;
			theta=Pij;
			Pij=Pji;
			Pji=theta;
			theta=Qij;
			Qij=Qji;
			Qji=theta;
			output<<setw(4)<<i<<setw(4)<<j<<setw(12)<<Pij<<setw(12)<<Qij<<setw(12)<<Pji<<setw(12)<<Qji<<endl;
		}
		else
		{
			output<<setw(4)<<i<<setw(4)<<j<<setw(12)<<Pij<<setw(12)<<Qij<<setw(12)<<Pji<<setw(12)<<Qji<<endl;
		}
	}
	output<<"有功损耗:"<<PLoss<<";无功损耗:"<<QLoss<<endl<<"系统最低电压:"<<Vmin<<";最低电压所在节点:"<<VminNode<<endl;
}






