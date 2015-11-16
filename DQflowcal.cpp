#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdio.h>
#include <iomanip.h>
#include <stdlib.h>
#include "definition.h"

#define Error 0.00001

////////////////////////
//////����ȫ�ֱ���//////
////////////////////////

/*****���ݶ�ȡ*****/
int					number;						//�û�����ֵ
char				file[10],*filename;			//�û������ļ��� �ļ���ָ��filename

/*****ϵͳ������Ϣ*****/
int					Mmax;						//����������
double				SB,V0;						//��׼����SB ��׼��ѹV0
int					N,Npv,Nb,Nt,Ng,Nl;			//�ڵ���N PV�ڵ���Npv ��·֧·��Nb ��ѹ��֧·��Nt ������ڵ���Ng ���ɽڵ���Nl
int					Nbg=0;						//�ӵ�֧·��
Branch_Info			*Branch;					//֧·����
Generator_Info		*Generator;					//������ڵ�����
Load_Info			*Load;						//���ɽڵ�����
Node_Info			*Node;						//�ڵ�����
PVNode_Info			*PVNode;					//PV�ڵ�����
NodalPower_Info		*NodalPower;				//�ڵ㹦��
NodalPower_Info		*GenePower;					//������ڵ㹦��
NodalVoltage_Info	*NodalVoltage;				//�ڵ��ѹ

/*****�ڵ��Ż�*****/
int					*OptNodeNum;				//�Ż��ڵ���������ָ��
int					ErrorNode;					//��������ڵ��

/*****���ɾ������*****/
int					*NYseq,*NYsum;				//���зǶԽ�Ԫ���׵�ַ����ָ��NYseq ���зǶԽ�Ԫ�ظ�������ָ��NYsum
U_Info				*U1,*U2;					//B'�����Ǿ���Ԫ��U1 B''�����Ǿ���Ԫ��U2
double				*D1,*D2;					//B'�ԽǾ���Ԫ��D1 B''�ԽǾ���Ԫ��D2
int					*NUsum1,*NUsum2;			//B'���ӱ������Ǿ�����зǶԽ�Ԫ���� B''���ӱ������Ǿ�����зǶԽ�Ԫ����
Yii_Info			*Yii;
Yij_Info			*Yij;						//���ɾ���Y�Խ�Ԫ��Yii �ǶԽ�Ԫ��Yij
Yii_Info			*Yii1;
Yij_Info			*Yij1;						//�����γ�B'�ľ���Y1�Խ�Ԫ��Yii1 �ǶԽ�Ԫ��Yij1
Yii_Info			*Yii2;
Yij_Info			*Yij2;						//�����γ�B''�ľ���Y2�Խ�Ԫ��Yii2 �ǶԽ�Ԫ��Yij2

double				NowError;
int					flag;
double				*DI1,*DI2;
double				Dtheta,Dv;
double				MaxError;
int					MaxErrorNode;
int					cnt;

////////////////////////
//////////����//////////
////////////////////////

/*****������*****/
void main()
{
	while(1)
	{
		/*step1����ʾ��ӭ��Ϣ*/

		cout<<endl<<"****************************PQ��������������****************************"<<endl<<endl;
		cout<<"			1.IEEE 5�ڵ�ϵͳ"<<endl;
		cout<<"			2.IEEE 14�ڵ�ϵͳ"<<endl;
		cout<<"			3.IEEE 30�ڵ�ϵͳ"<<endl;
		cout<<"			4.IEEE 57�ڵ�ϵͳ"<<endl;
		cout<<"			5.IEEE 118�ڵ�ϵͳ"<<endl;
		cout<<"			6.�Զ����ļ�"<<endl;
		cout<<"��������Ų��س���";



		/*step2���û�����*/

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
			cout<<"�������ļ�����";
			cin>>file;
			filename=file;
			break;
		}



		/*step3��ȷ���ļ����Ƿ���ȷ ���ļ�����ȡ����*/

		ifstream in(filename,ios::nocreate);
		if(!in)
		{
			cout<<"�ļ������󣡰����������"<<endl;
			system("pause");
			continue;
		}
		else
			Read(filename);



		/*step3��T2�ڵ����Ż�*/

		NodeOpt();



		/*step4���γɵ��ɾ���*/

		YFormation();



		/*step5���γ�B' B''�����ӱ�*/

		NUsum1=new int[N+1];						
		NUsum2=new int[N+1];
		D1=new double[N+1];
		D2=new double[N+1];
		U1=new U_Info[(Nb+Nt)*8];
		U2=new U_Info[(Nb+Nt)*8];
		BFormation(Yii1,Yij1,1,NUsum1,U1,D1);
		BFormation(Yii2,Yij2,2,NUsum2,U2,D2);

		//���
		cout<<endl<<"5.�������ӱ�"<<endl;
		//��������Ǿ�����зǶԽ�Ԫ����
		cout<<endl<<"��һ���ӱ������Ǿ�����зǶԽ�Ԫ������"<<endl;
		for(int i=1;i<N;i++)
			cout<<"��"<<i<<"��	"<<NUsum1[i]<<"��"<<endl;
		cout<<endl<<"�ڶ����ӱ������Ǿ�����зǶԽ�Ԫ������"<<endl;
		for(i=1;i<N;i++)
			cout<<"��"<<i<<"��	"<<NUsum2[i]<<"��"<<endl;
		cout<<endl;
		system("pause");
		//������ӱ�Խ�Ԫ��
		cout<<endl<<"���ӱ�Խ�Ԫ��D1:"<<endl;
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
		cout<<endl<<"���ӱ�Խ�Ԫ��D2:"<<endl;
		for(i=1;i<N;i++)
			cout<<i<<"	"<<D2[i]<<endl;
		cout<<endl;
		system("pause");
		//������ӱ�������Ԫ
		cout<<endl<<"��һ���ӱ�������ԪU1="<<endl;
		for(i=1;i<=t1;i++)
			cout<<i<<setw(12)<<U1[i].value<<endl;
		cout<<endl<<"�ڶ����ӱ�������ԪU2="<<endl;
		for(i=1;i<=t2;i++)
			cout<<i<<setw(12)<<U2[i].value<<endl;
		cout<<endl;
		system("pause");



		/*step6�����ڵ��ѹ��ֵ*/

		NodalVoltage=new NodalVoltage_Info[N+1];
		Initial(1,0);



		/*step7�����е���*/
		cnt=1;
		NodalPower=new NodalPower_Info[N+1];
		GenePower=new NodalPower_Info[N+1];
		DI1=new double[N+1];
		DI2=new double[N+1];
		NowError=100;									//������ֵ
		cout<<endl<<"���������ʼ����"<<endl<<endl;
		system("pause");
		while(NowError>Error)							//������Ҫ�����������
		{
			//P����
			PQCal(1);									//����ڵ㹦��
			PQErrorCal(DI1,1);							//����ڵ㹦�����
			SolvEqu(U1,D1,DI1,NUsum1);					//�����Է�����
			for(i=1;i<N;i++)							//����
			{
				Dtheta=DI1[i]/NodalVoltage[i].V;
				NodalVoltage[i].theta-=Dtheta;
			}
			NowError=MaxError;
			MaxErrorNode=ErrorNode;

			//Q����
			PQCal(2);									//����ڵ㹦��
			PQErrorCal(DI2,2);							//����ڵ㹦�����
			SolvEqu(U2,D2,DI2,NUsum2);					//�����Է�����
			for(i=1;i<N;i++)							//����
			{
				Dv=DI2[i];
				NodalVoltage[i].V-=Dv;
			}

			if(NowError<MaxError)						//�ж�������
			{
				NowError=MaxError;
				MaxErrorNode=ErrorNode;
			}

			//���ÿ�ε������
			cout<<"============================================================="<<endl;
			cout<<endl<<"��"<<cnt<<"�ε���"<<endl;
			cout<<"�����"<<NowError<<"	������ڵ�ţ�"<<MaxErrorNode<<endl<<endl;
			for(i=1;i<N;i++)
				cout<<"�ڵ�"<<i<<"	��P:"<<setw(15)<<DI1[i]<<setw(5)<<"��:"<<setw(15)<<NodalVoltage[i].theta<<endl;
			for(i=1;i<N;i++)
				cout<<"�ڵ�"<<i<<"	��Q:"<<setw(15)<<DI2[i]<<setw(5)<<"V:"<<setw(15)<<NodalVoltage[i].V<<endl;
			cnt++;
			cout<<endl;
			if(NowError<Error)
			{
				cout<<"��"<<cnt<<"�ε���������"<<endl<<endl;
			}
			if(cnt>Mmax)
			{
				cout<<"����������������!"<<endl;
				break;
			}
		}



		/*step8��������*/
		OutData();
		system("pause");
	}
}



/*****���ݶ�ȡ*****/
void Read(char *filename)
{
	char head[100],enter;						//��ͷ ���з��Ĵ���



	/*step1�����ļ�*/

	ifstream in(filename);



	/*step2����ȡ��Ϣ*/

	//������Ϣ
	in.getline(head,100);						//������ͷ
	in>>Mmax>>SB>>V0;							//��ȡ����������Mmax ϵͳ��׼����SB ϵͳ��׼��ѹV0
	in>>enter;									//�������з�
	in.getline(head,100);						//������ͷ
	in>>N>>Nb>>Nt>>Ng>>Nl;						//��ȡ�ڵ���N ��·֧·��Nb ��ѹ��֧·��Nt ������ڵ���Ng ���ɽڵ���Nl

	//֧·
	in>>enter;									//�������з�
	in.getline(head,100);						//������ͷ
	Branch=new Branch_Info[Nb+Nt+1];
	for(int i=1;i<=Nb+Nt;i++)
	{
		in>>Branch[i].i>>Branch[i].j>>Branch[i].R>>Branch[i].X>>Branch[i].YK;
	}

	//������ڵ�
	in>>enter;									//�������з�
	in.getline(head,100);						//������ͷ
	Generator=new Generator_Info[Ng+1];
	for(i=1;i<=Ng;i++)
	{
		in>>Generator[i].i>>Generator[i].j>>Generator[i].P>>Generator[i].Q>>Generator[i].V;
	}

	//���ɽڵ�
	Load=new Load_Info[Nl+1];
	for(i=1;i<=Nl;i++)
	{
		in>>Load[i].i>>Load[i].j>>Load[i].P>>Load[i].Q>>Load[i].V;
	}
    in.close();



	/*step3:���*/

	cout<<endl<<"***************************************************************"<<endl;
	cout<<endl<<"1.ϵͳ��Ϣ"<<endl;

	//������Ϣ
	cout<<endl<<"=======================������Ϣ======================="<<endl;
	cout<<"����������\t"<<"ϵͳ��׼����\t"<<"ϵͳ��׼��ѹ"<<endl;
	cout<<Mmax<<"\t\t"<<SB<<"\t\t"<<V0<<endl;
	cout<<"ϵͳ�ڵ���\t"<<"��·֧·��\t"<<"��ѹ��֧·��\t"<<"������ڵ���\t"<<"���ɽڵ���"<<endl;
	cout<<N<<"\t\t"<<Nb<<"\t\t"<<Nt<<"\t\t"<<Ng<<"\t\t"<<Nl<<endl<<endl;
	system("pause");

	//֧·��Ϣ
	cout<<endl<<"=======================֧·��Ϣ======================="<<endl;
	cout<<"֧·���ӵ�i\t"<<"֧·���ӵ�j\t"<<"֧·����\t"<<"֧·�翹\t"<<"֧·�ӵص���"<<endl;
	for(i=1;i<=Nb+Nt;i++)
	{
		cout<<Branch[i].i<<"\t\t"<<Branch[i].j<<"\t\t"<<Branch[i].R<<"\t\t"<<Branch[i].X<<"\t\t"<<Branch[i].YK<<endl;
	}
	cout<<endl;
	system("pause");

	/*�ڵ���Ϣ*/
	cout<<endl<<"=======================�ڵ���Ϣ======================="<<endl;
	cout<<"�ڵ���\t"<<"�ڵ�����\t"<<"�ڵ��й�����\t"<<"�ڵ��޹�����\t"<<"�ڵ��ѹ"<<endl;
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

/*****�ڵ����Ż�*****/
void NodeOpt()
{
	Node=new Node_Info[N+1];					//�ڵ�ṹ��
	int Balance=0;								//ƽ��ڵ��



	/*step1����ʼ��*/

	for(int i=1;i<=N;i++)
	{
		Node[i].num=i;
		Node[i].outnum=0;
		for(int ii=1;ii<=20;ii++)
			Node[i].connect[ii]=0;
	}



	/*step2��ȷ��ƽ��ڵ� �ҽ����������Ϊ1000*/

	for(i=1;i<=Ng;i++)
	{
		if(Generator[i].j==0)
			Balance=Generator[i].i;
	}
	Node[Balance].outnum=1000;



	/*step3��ȷ�����ڵ����ӹ�ϵ*/

	for(i=1;i<=Nb+Nt;i++)            
	{
		//�ж��Ƿ�Ϊ�ӵ�֧· ������������ӹ�ϵ
		if(Branch[i].i==Branch[i].j)
		{
			Nbg++;
			continue;
		}

		//�ж��Ƿ�Ϊ����ƽ��ڵ��֧· ������������ӹ�ϵ
		else if(Branch[i].i==Balance||Branch[i].i==Balance)		
			continue;

		//�ж��Ƿ��������ӹ�ϵ ��������������ӹ�ϵ
		for(int ii=1;ii<=Node[Branch[i].i].outnum;ii++)
		{
			if(Node[Branch[i].i].connect[ii]==Branch[i].j)
				break;
		}

		//������ӹ�ϵ
		if(ii>Node[Branch[i].i].outnum)							
		{
			Node[Branch[i].i].outnum++;
			Node[Branch[i].i].connect[Node[Branch[i].i].outnum]=Branch[i].j;
			Node[Branch[i].j].outnum++;
			Node[Branch[i].j].connect[Node[Branch[i].j].outnum]=Branch[i].i;
		}
	}



	/*step4������PV�ڵ���Ϣ*/

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



	/*step5��������ڵ��������*/

	cout<<endl<<"***************************************************************"<<endl;
	cout<<endl<<"2.�ڵ��������"<<endl<<endl;

	for(i=1;i<=N;i++)
	{
		if(Node[i].outnum>0&&i!=Balance)
		{
			cout<<"�ڵ�"<<i<<"������Ϊ"<<Node[i].outnum<<"�����ӵĽڵ��У�";
			for(int ii=1;ii<=Node[i].outnum;ii++)
				cout<<Node[i].connect[ii]<<" ";
			cout<<endl;
		}
	}
	cout<<endl<<"PV�ڵ���"<<Npv<<"����Ϊ";
	for(i=1;i<=Npv;i++)
	{
		cout<<PVNode[i].i<<" ";
	}
	cout<<endl<<"ƽ��ڵ�Ϊ"<<Balance<<endl;
	cout<<"�ӵ�֧·��"<<Nbg<<"��"<<endl<<endl;
	system("pause");



	/*step6���붯̬�ڵ����Ż�*/

	OptNodeNum=new int[N+1];
	int leastnum;
	OptNodeNum[N]=Balance;
	for(i=1;i<=N-1;i++)
	{
		//�ҳ���������С�ڵ�
		leastnum=1;
		for(int j=1;j<=N;j++)								
		{
			if(Node[j].outnum<Node[leastnum].outnum)
				leastnum=j;
		}
		OptNodeNum[i]=leastnum;												//���½ڵ����д���Ż��������
		/*
		cout<<"�˴���ȥ�ڵ�Ϊ��"<<leastnum<<endl;							//�����ȥ�ڵ�
		system("pause");
		*/

		//��ȥ�ڵ� ���³�����
		for(j=1;j<=Node[leastnum].outnum;j++)
		{
			for(int k=1;k<=Node[Node[leastnum].connect[j]].outnum;k++)
				if(Node[Node[leastnum].connect[j]].connect[k]==leastnum)	//�ж���ȥ�ڵ��Ƿ��Ǹýڵ������б�����һ��
					break;
			if(k==Node[Node[leastnum].connect[j]].outnum)					//����ȥ�ڵ��Ǹýڵ������б�����һ�� ��ӽ��ýڵ������-1
				Node[Node[leastnum].connect[j]].outnum--;
			else															//����ȥ�ڵ㲻�Ǹýڵ������б�����һ�� �򽫸ýڵ������-1 �������������ڵ�ǰ��
			{
			for(k=k;k<Node[Node[leastnum].connect[j]].outnum;k++)
				Node[Node[leastnum].connect[j]].connect[k]=Node[Node[leastnum].connect[j]].connect[k+1];
			Node[Node[leastnum].connect[j]].outnum--;
			}
		}
		/*
		for(j=1;j<=N;j++)												//�����ȥ��ڵ����ӱ仯���
		{
			if(Node[j].outnum!=1000&&j!=leastnum)
			{
				cout<<"�˴���ȥ��ڵ�"<<j<<"�ĳ�����Ϊ"<<Node[j].outnum<<"���������ӵĽڵ��У�";
				for(int k=1;k<=Node[j].outnum;k++)
				cout<<"�ڵ�"<<Node[j].connect[k]<<" ";
			cout<<endl;
			}
		}
		system("pause");
		*/

		//׷��֧· ���³�����
		for(j=1;j<=Node[leastnum].outnum-1;j++)							//�ж�����ȥ�ڵ�leastnum�����ĵ�j�ڵ�͵�k�ڵ��Ƿ����� ������׷��֧·
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
		for(j=1;j<=N;j++)												//���׷�Ӻ�ڵ����ӱ仯���
		{
			if(Node[j].outnum!=1000&&j!=leastnum)
			{
				cout<<"�˴�׷�Ӻ�ڵ�"<<j<<"�ĳ�����Ϊ"<<Node[j].outnum<<"���������ӵĽڵ��У�";
				for(int k=1;k<=Node[j].outnum;k++)
					cout<<"�ڵ�"<<Node[j].connect[k]<<" ";
				cout<<endl;
			}
		}
		system("pause");
		*/
		Node[leastnum].outnum=1000;										//����������ڵ��ų��Ƚ�
	}



	/*step7�����½ڵ���*/

	for(i=1;i<=N;i++)
		Node[OptNodeNum[i]].num=i;

	//PV�ڵ�
	for(i=1;i<=Npv;i++)
		PVNode[i].i=Node[PVNode[i].i].num;
	//������ڵ�
	for(i=1;i<=Ng;i++)
		Generator[i].i=Node[Generator[i].i].num;
	//���ɽڵ�
	for(i=1;i<=Nl;i++)
		Load[i].i=Node[Load[i].i].num;


	/*step8���ڵ���������*/
	//PV�ڵ�
	PVNode_Info PVtemp;
	for(i=1;i<Npv;i++)
		for(int j=Npv;j>i;j--)
			if(PVNode[j-1].i>PVNode[j].i)
			{
				PVtemp=PVNode[j];
				PVNode[j]=PVNode[j-1];
				PVNode[j-1]=PVtemp;
			}
	//������ڵ�
	Generator_Info Genetemp;
	for(i=1;i<Ng;i++)
		for(int j=Ng;j>i;j--)
			if(Generator[j-1].i>Generator[j].i)
			{
				Genetemp=Generator[j];
				Generator[j]=Generator[j-1];
				Generator[j-1]=Genetemp;
			}
	//���ɽڵ�
	Load_Info Loadtemp;
	for(i=1;i<Nl;i++)
		for(int j=Nl;j>i;j--)
			if(Load[j-1].i>Load[j].i)
			{
				Loadtemp=Load[j];
				Load[j]=Load[j-1];
				Load[j-1]=Loadtemp;
			}



	/*step9������֧·�ͽڵ���Ϣ*/

	//��·֧·
	for(i=1;i<=Nb;i++)													
	{
		Branch[i].i=Node[Branch[i].i].num;
		Branch[i].j=Node[Branch[i].j].num;
	}
	//��ѹ��֧·
	for(i=Nb+1;i<=Nb+Nt;i++)
	{
		Branch[i].i=Node[Branch[i].i].num;
		Branch[i].j=-Node[Branch[i].j].num;								//������Ϊ��ѹ��֧·��־
	}
	//֧·��Ÿ���
	Branch_Info Branchtem;												//��ʱ����
	int k=0;
	for(i=1;i<=Nb+Nt-Nbg;i++)         
	{
		if(abs(Branch[i].i)>abs(Branch[i].j))							//֧·���˽ڵ�Ž�С������ǰ��
		{
			int j=Branch[i].i;
			Branch[i].i=Branch[i].j;
			Branch[i].j=j;

		}
		else if(Branch[i].i==Branch[i].j)								//���ӵ�֧·�������
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
	for(i=1;i<=N;i++)													//��֧·���սڵ��Ŵ�С��������
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



	/*step10������Ż����*/

	cout<<endl<<"***************************************************************"<<endl;
	cout<<endl<<"3.T2�ڵ����Ż�"<<endl<<endl;

	//�ڵ��Ż����
	for(i=1;i<=N;i++)
		cout<<"ԭ��"<<i<<"�ڵ㣬���붯̬�Ż�����Ϊ"<<Node[i].num<<endl;
	cout<<endl;
	system("pause");
	//�Ż���֧·��Ϣ
	cout<<endl<<"=======================�Ż���֧·��Ϣ======================="<<endl;
	cout<<"֧·���ӵ�i\t"<<"֧·���ӵ�j\t"<<"֧·����\t"<<"֧·�翹\t"<<"֧·�ӵص���"<<endl;
	for(i=1;i<=Nb+Nt;i++)
	{
		cout<<Branch[i].i<<"\t\t"<<Branch[i].j<<"\t\t"<<Branch[i].R<<"\t\t"<<Branch[i].X<<"\t\t"<<Branch[i].YK<<endl;
	}
	cout<<endl;
	system("pause");
	//�Ż���ڵ���Ϣ
	cout<<endl<<"=======================�Ż���ڵ���Ϣ======================="<<endl;
	cout<<"�ڵ���\t"<<"�ڵ�����\t"<<"�ڵ��й�����\t"<<"�ڵ��޹�����\t"<<"�ڵ��ѹ"<<endl;
	for(i=1;i<=Ng;i++)
	{
		cout<<Generator[i].i<<"\t\t"<<Generator[i].j<<"\t\t"<<Generator[i].P<<"\t\t"<<Generator[i].Q<<"\t\t"<<Generator[i].V<<endl;
	}
	for(i=1;i<=Nl;i++)
	{
		cout<<Load[i].i<<"\t\t"<<Load[i].j<<"\t\t"<<Load[i].P<<"\t\t"<<Load[i].Q<<"\t\t"<<Load[i].V<<endl;
	}
	cout<<endl;
	cout<<"PV�ڵ���"<<Npv<<"����Ϊ";
	for(i=1;i<=Npv;i++)
	{
		cout<<PVNode[i].i<<" ";
	}
	cout<<endl;
	cout<<"ƽ��ڵ�Ϊ�ڵ�"<<N<<endl;
	cout<<"�ӵ�֧·��"<<Nbg<<"��"<<endl<<endl;
	system("pause");
}

/*****���ɾ����γ�*****/
void YFormation()
{
	int i,j;															//���ӽڵ���
	double R,X,YK,Zmag2;												//֧·���� �翹 �ӵص��� �迹��ֵ
	double G,B,b;														//����ʵ�� �鲿
	Yii=new Yii_Info[N+1];
	Yij=new Yij_Info[Nb+Nt+1];											//Y��ŵ��ɾ���
	Yii1=new Yii_Info[N+1];
	Yij1=new Yij_Info[Nb+Nt+1];											//Y1�����γ�B'�������ǽӵ�֧·��
	Yii2=new Yii_Info[N+1];
	Yij2=new Yij_Info[Nb+Nt+1];											//Y2�����γ�B''����������·���裩
	NYsum=new int[N+1];
	NYseq=new int[N+1];



	/*step1����ʼ��*/
	for(int n=1;n<=N;n++)												//��ʼ���Խ�Ԫ��
	{
		Yii[n].G=0;
		Yii[n].B=0;
		Yii1[n].G=0;
		Yii1[n].B=0;
		Yii2[n].G=0;
		Yii2[n].B=0;
		NYsum[n]=0;
	}
	for(n=1;n<=Nb+Nt;n++)											//��ʼ���ǶԽ�Ԫ��
	{
		Yij[n].G=0;
		Yij[n].B=0;
		Yij1[n].G=0;
		Yij1[n].B=0;
		Yij2[n].G=0;
		Yij2[n].B=0;
	}



	/*step2���γɲ����ǽӵ�֧·�ĵ��ɾ���*/
	for(n=1;n<=Nb+Nt;n++)
	{
		i=abs(Branch[n].i);
		j=abs(Branch[n].j);
		R=Branch[n].R;
		X=Branch[n].X;
		YK=Branch[n].YK;

		if(i==j)														//�����Ե�֧·
			continue;

		//�����֧·����
		Zmag2=R*R+X*X;													//�迹��ֵ
		G=R/Zmag2;														//����ʵ��
		B=-X/Zmag2;														//�����鲿
		b=-1/X;															//�����ǵ���ĵ����鲿

		if(Branch[n].i<0||Branch[n].j<0)
		{
			Yij[n].G=-G/YK;												//�����ѹ��֧·������ ���ǶԽ�Ԫ��
			Yij[n].B=-B/YK;
			Yij1[n].G=-G/YK;
			Yij1[n].B=-B/YK;
			Yij2[n].G=0;
			Yij2[n].B=-b/YK;

			Yii[i].G=Yii[i].G+G/YK;										//�����ѹ���ڵ��Ե��� ���Խ�Ԫ��
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
			Yij[n].G=-G;												//������·֧·������ ���ǶԽ�Ԫ��
			Yij[n].B=-B;
			Yij1[n].G=-G;
			Yij1[n].B=-B;
			Yij2[n].G=0;
			Yij2[n].B=-b;

			Yii[i].G=Yii[i].G+G;										//������·�ڵ��Ե��� ���Խ�Ԫ��
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



	/*step3��׷�ӽӵ�֧·*/
	for(n=1;n<=Nb+Nt;n++)
	{
		i=Branch[n].i;
		j=Branch[n].j;
		R=Branch[n].R;
		X=Branch[n].X;
		YK=Branch[n].YK;

		if(i==j)													//����Ե�֧·
		{
			Yii[i].B=Yii[i].B+1.0/X;
			Yii2[i].B=Yii2[i].B+1.0/X;
			continue;
		}

		if(i<0||j<0)
		{
			if(i<0)
			{
				i=abs(i);											//iΪ�Ǳ�׼��Ȳ�
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
				j=abs(j);											//jΪ�Ǳ�׼��Ȳ�
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
			B=YK;													//��·
			b=YK;
			Yii[i].B=Yii[i].B+B;
			Yii[j].B=Yii[j].B+B;
			Yii2[i].B=Yii2[i].B+b;
			Yii2[j].B=Yii2[j].B+b;
		}
	}



	/*step4��������ɾ���*/
	cout<<endl<<"4.���㵼�ɾ���"<<endl;
	//���Y
	cout<<endl<<"ϡ�赼�ɾ���Խ�Ԫ:Yii="<<endl;
	for(i=1;i<=N;i++)
		cout<<i<<setw(10)<<Yii[i].G<<"+"<<setw(10)<<Yii[i].B<<"j"<<endl;
	cout<<"ϡ�赼�ɾ���ǶԽ�Ԫ:Yij="<<endl;
	for(i=1;i<=Nb+Nt-Nbg;i++)
	{
		cout<<i<<setw(10)<<Yij[i].G<<"+"<<setw(10)<<Yij[i].B<<"j"<<endl;
	}
	cout<<endl;
	system("pause");
	//���Y1
	cout<<endl<<"ϡ�赼�ɾ���Խ�Ԫ:Yii1="<<endl;
	for(i=1;i<=N;i++)
		cout<<i<<setw(10)<<Yii1[i].G<<"+"<<setw(10)<<Yii1[i].B<<"j"<<endl;
	cout<<"ϡ�赼�ɾ���ǶԽ�Ԫ:Yij1="<<endl;
	for(i=1;i<=Nb+Nt-Nbg;i++)
	{
		cout<<i<<setw(10)<<Yij1[i].G<<"+"<<setw(10)<<Yij1[i].B<<"j"<<endl;
	}
	cout<<endl;
	system("pause");
	//���Y2
	cout<<endl<<"ϡ�赼�ɾ���Խ�Ԫ:Yii2="<<endl;
	for(i=1;i<=N;i++)
		cout<<i<<setw(10)<<Yii2[i].G<<"+"<<setw(10)<<Yii2[i].B<<"j"<<endl;
	cout<<"ϡ�赼�ɾ���ǶԽ�Ԫ:Yij="<<endl;
	for(i=1;i<=Nb+Nt-Nbg;i++)
	{
		cout<<i<<setw(10)<<Yij2[i].G<<"+"<<setw(10)<<Yij2[i].B<<"j"<<endl;
	}
	cout<<endl;
	system("pause");
}

/*****�γ�B' B''�����ӱ�*****/
void BFormation(Yii_Info *Yii,Yij_Info *Yij,int flag,int *NUsum,U_Info *U,double *D)
//����Խ�Ԫ��Yii �ǶԽ�Ԫ��Yij Yii1��Yij1�����γ�B' Yii2��Yij2�����γ�B''
//��ʶ����flag=1�γ�B'���ӱ� flag=2�γ�B''���ӱ�
//NUsum��������Ǿ�����зǶԽ�Ԫ����
//D������ӱ�ĶԽ�Ԫ��
{
	int n_pv,i_pv;													//�ڵ��������n_pv �ڵ��i_pv
	int i,j;														//�кű���i ���±����j
	int n,n_u;														//��ʱ��������n ���ӱ������Ǿ���Ԫ�ؼ�������
	int i_above;													//��ȥ�к�i��������
	double	*B;														//ϵ������BԪ��
	double Btemp;
	B=new double[N+1];												//��ʱ����
	n_pv=1;
	i_pv=PVNode[1].i;

	for(i=1;i<N;i++)
	{
		if(flag==2&&i==i_pv)
		{
			//F����
			n_pv++;
			i_pv=PVNode[n_pv].i;
			NUsum[i]=0;
			D[i]=0.0;
		}
		else
		{
			//A����
			for(n=i+1;n<N;n++)
				B[n]=0;
			B[i]=Yii[i].B;

			//B����
			for(n=NYseq[i];n<NYseq[i+1];n++)
			{
				j=Yij[n].j;
				B[j]=Yij[n].B;
			}

			//C����
			if(flag==2)
			{
				for(n=1;n<=Npv;n++)
					B[PVNode[n].i]=0;
			}

			//D����
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

			//E����
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



/*���ڵ��ѹ��ֵ*/
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



/*�ڵ㹦�ʼ���*/
void PQCal(int flag)							//flag=1����P-theta���� flag=2����Q-V����
{
	double Vi,A,B,vv,theta;						//i�ڵ��ѹ��ֵVi ��ʱ����A B vv theta
	int j;										//��ʱ����j
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
		if(i==N)								//���һ���޷ǶԽ�Ԫ������
			break;
		for(int n=NYseq[i];n<NYseq[i+1];n++)	//�ҵ���i�еķǶԽ�Ԫ��
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

/*�����������*/
void PQErrorCal(double *DI,int flag)								//�������DI flag=1����P��� flag=2����Q���
{
	MaxError=0;												//��������MaxError
	double Vi;														//i���ѹ��ֵVi
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

		//����
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

		//�ڵ�ע�빦��
		if(flag==1)
			Wi=Wi-NodalPower[i].P;
		else
			Wi=Wi-NodalPower[i].Q;

		//�����
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

		if(i==N)									//ƽ��ڵ㲻���㹦�����
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

/*�����Է�����*/
void SolvEqu(U_Info *U,double *D,double *DI,int *NUsum)
{
	int	n_u=1;
	int i,j,count;
	double	DItemp;
	for(i=1;i<N;i++)
	{
		DItemp=DI[i];
		for(count=1;count<=NUsum[i];count++)				//ǰ��
		{
			j=U[n_u].j;
			DI[j]=DI[j]-DItemp*U[n_u].value;
			n_u++;
		}
		DI[i]=DItemp*D[i];
	}
	for(i=N-1;i>0;i--)										//�ش�
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

///////////////////////////��������Ӻ���////////////////////////////////////////
void OutData()
{
	ofstream	output;											//���
	int		i,j,k,n,i_g,n_g,VminNode;
	double	theta,Vmin,V,P,Q;
	double	PLoss,QLoss,R,X,YK,Vi,Vj,Ei,Ej,Fi,Fj,DE,DF,Ir,Ii,Pij,Qij,Pji,Qji,Zmag;
	Vmin=NodalVoltage[1].V;           //ϵͳ��͵�ѹ      
	i_g=Generator[1].i;       
	VminNode=1;					//ϵͳ��͵�ѹ ���ڽڵ��                     
	n_g=1;
	output.open("output.txt");
	///////////////////////////���///////////////////////////////
	output<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(6); //��Ļ�����ʽ���ơ�
	output<<"******************************�ڵ���Ϣ�����*******************************"<<endl;
	output<<setw(4)<<"�ڵ��"<<setw(12)<<"�ڵ��ѹV"<<setw(12)<<"��Ǧ�"<<setw(12)<<"������й�"<<setw(12)<<"������޹�"<<setw(12)<<"PL"<<setw(12)<<"QL"<<endl;
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
	output<<endl<<"**********************֧·���������**************************"<<endl;
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
	output<<"�й����:"<<PLoss<<";�޹����:"<<QLoss<<endl<<"ϵͳ��͵�ѹ:"<<Vmin<<";��͵�ѹ���ڽڵ�:"<<VminNode<<endl;
}






