
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"


// ????? Why trouble use const int** and const double**

// ##### =>
#define RouteLimit 11 //The first one is length. RouteLimit = Routelength + 1.
#define CabsNum 96
#define AccRecordPeriod 1000
// <= #####

int AccNum;

// Function Series 1: Basic

// Copy a to b.
void CopyInt(const int* a, int* b, const int n){
    int i;
    
    for(i=0;i<n;i++)
        b[i]=a[i];
}

void SwapInt(int* a, int* b){
    int temp;
    
    temp=*a;
    *a=*b;
    *b=temp;
}

void SortInt_NoPlace(int* Value, const int size){
    int i,j;
    
    for(i=0;i<size-1;i++){
        for(j=0;j<size-1-i;j++){
            if(Value[j]>Value[j+1]){
                SwapInt(&Value[j],&Value[j+1]);
            }
        }
    }
}

void Route2RouteVec(int** Route, int* RouteVec, const int CabsNumPerCore){
    int i,j,k;
    
    for(i=0,k=0;i<CabsNumPerCore;i++)
        for(j=0;j<=Route[i][0];j++)
            RouteVec[k++]=Route[i][j];
}

void RouteVec2Route(int** Route, const int* RouteVec, const int CabsNumPerCore){
    int i,j,k;
    
    for(i=0,k=0;i<CabsNumPerCore;i++){
        Route[i][0]=RouteVec[k++];
        for(j=1;j<=Route[i][0];j++)
            Route[i][j]=RouteVec[k++];
    }
}

double TwoVec_SimilarityRate(const int* a, const int* b, const int n){
    int i,SimNum;
    
    for(i=0,SimNum=0;i<n;i++){
        if(a[i]==b[i])
            SimNum++;
    }
    return ((double)SimNum)/(double)n;
}

struct pair_ERank{
    double E;
    int rank;
};

// Function Series 2: SA OneX

double PTD_OneX(const int* Route1, const double DInf, double** SpotsPos, const double* PosCab, double* DistVec, double* RouteProb){
    int i,j;
    double CumuRouteLength,CumuProb;
    double PTD_E;
    
    DistVec[0]=pow(pow(PosCab[0]-SpotsPos[Route1[1]][0],2)+pow(PosCab[1]-SpotsPos[Route1[1]][1],2),0.5);
    for(i=1;i<Route1[0];i++){
        DistVec[i]=pow(pow(SpotsPos[Route1[i]][0]-SpotsPos[Route1[i+1]][0],2)+pow(SpotsPos[Route1[i]][1]-SpotsPos[Route1[i+1]][1],2),0.5);
    }
    for(i=0;i<Route1[0];i++){
        RouteProb[i]=SpotsPos[Route1[i+1]][2];
    }
    PTD_E=0;
    CumuRouteLength=0;
    CumuProb=1;
    CumuRouteLength+=DistVec[0];
    PTD_E+=CumuRouteLength*RouteProb[0];
    CumuProb*=1-RouteProb[0];
    
    for(i=1;i<Route1[0];i++){
        CumuRouteLength+=DistVec[i];
        PTD_E+=CumuRouteLength*CumuProb*RouteProb[i];
        CumuProb*=1-RouteProb[i];
    }
    PTD_E+=DInf*CumuProb;//Here the DInf actually incorporates D1+D2+... .
    return PTD_E;
}

// 1 is accept and 0 is not.
int SA_AcceptOrNot(const double PTD_E0, const double PTD_E1, const double Temperature){
    int AcceptFlag;
    double AccProb;
    
    if(PTD_E1<PTD_E0)
        AcceptFlag=1;
    else{
        if(Temperature<0.000000000001)
            AcceptFlag=0;
        else if((PTD_E1-PTD_E0)/Temperature>30)
            AcceptFlag=0;
        else{
            AccProb=exp((PTD_E0-PTD_E1)/Temperature);
            if(rand()/RAND_MAX<AccProb)
                AcceptFlag=1;
            else
                AcceptFlag=0;
        }
    }
    return AcceptFlag;
}

// Function Series 3: SA MulX

double PTD_MulX(int** MulX_Route1, const double* MulX_DInf, double** SpotsPos, double** MulX_PosCab, const int XNum, double* DistVec, double* RouteProb){
    int i;
    double PTD_sum;
    
    PTD_sum=0.0;
    for(i=0;i<XNum;i++){
        PTD_sum+=PTD_OneX(MulX_Route1[i],MulX_DInf[i],SpotsPos,MulX_PosCab[i],DistVec,RouteProb);
    }
    return PTD_sum;
}

void SA_MulXIni_NoCandiSeq(int** SoluMulX, const int XNum, const int SpotsNum){
    int i,j,k,temp;
    int *Seq;
    Seq=(int*)malloc(SpotsNum*sizeof(int));
    
    for(i=0;i<SpotsNum;i++)
        Seq[i]=i;
    for(i=0,k=0;i<XNum;i++){
        for(j=1;j<=SoluMulX[i][0];j++){
            temp=rand()%(SpotsNum-k);
            SwapInt(&Seq[k],&Seq[k+temp]);
            SoluMulX[i][j]=Seq[k];
            k++;
        }
    }
    
    free(Seq);
}

void SA_GloMove_MulXSwap_NoCandiSeq(int** SoluMulX_new, const int XNum, const int SpotsNum, int* rr){
    int i,j,i0,j0,flag,temp;
    
    flag=0;
    i0=rand()%XNum;
    j0=1+rand()%SoluMulX_new[0][0];// All route length should be the same.
    temp=rand()%SpotsNum;
    for(i=0;i<XNum;i++){
        for(j=1;j<=SoluMulX_new[i][0];j++){
            if(SoluMulX_new[i][j]==temp){
                flag=1;
                break;
            }
        }
        if(flag==1)
            break;
    }
    if(flag){
        SwapInt(&SoluMulX_new[i0][j0],&SoluMulX_new[i][j]);
        rr[1]=i;
    }else{
        SoluMulX_new[i0][j0]=temp;
        rr[1]=-1;
    }
    rr[0]=i0;
}

void SA_GloMove_MulXSwap_PMode2(int** SoluMulX_new, const int XNum, int* CandiSeq, const int CandiNum, const int* OriSoluVec, int OriSoluLength, int** SoluRecord, int* rr){
    int i,j,flag,i0,j0,RandCandi;
    
    flag=0;
    i0=rand()%XNum;
    j0=1+rand()%SoluMulX_new[0][0];// All route length should be the same.
    RandCandi=rand()%(OriSoluLength+CandiNum);
    for(i=0;i<XNum;i++){
        for(j=1;j<=SoluMulX_new[i][0];j++){
            if(SoluRecord[i][j]==RandCandi){
                flag=1;
                break;
            }
        }
        if(flag==1)
            break;
    }
    if(flag){
        SwapInt(&SoluRecord[i0][j0],&SoluRecord[i][j]);
        SwapInt(&SoluMulX_new[i0][j0],&SoluMulX_new[i][j]);
        rr[1]=i;
    }else{
        SoluRecord[i0][j0]=RandCandi;
        if(RandCandi<OriSoluLength)
            SoluMulX_new[i0][j0]=OriSoluVec[RandCandi];
        else{
            SoluMulX_new[i0][j0]=CandiSeq[RandCandi-OriSoluLength];
        }
        rr[1]=-1;
    }
    rr[0]=i0;
}

// Each core has its own CandiSeq. It needs to initialize its own SoluMulX and shorten the CandiNum.
// !!!!! The memory of the tail of CandiSeq is wasted.
void SA_MulXIni(int** SoluMulX, const int XNum, int* CandiSeq, int* CandiNumAdd){
    int i,j,k,temp;
    
    for(i=0,k=0;i<XNum;i++){
        for(j=1;j<=SoluMulX[i][0];j++){
            temp=rand()%(*CandiNumAdd-k);
            SwapInt(&CandiSeq[*CandiNumAdd-1-k],&CandiSeq[temp]);
            SoluMulX[i][j]=CandiSeq[*CandiNumAdd-1-k];
            k++;
        }
    }
    *CandiNumAdd=*CandiNumAdd-k;
}

void SA_GloMove_MulXSwap(int** SoluMulX_new, const int XNum, int* CandiSeq, const int CandiNum, int* r){
    int i,SoluMulX_num,RandIndex1,RandIndex2,a10,a11,b10,b11;
    
    SoluMulX_num=0;
    for(i=0;i<XNum;i++){
        SoluMulX_num+=SoluMulX_new[i][0];
    }
    RandIndex1=rand()%SoluMulX_num;
    RandIndex2=rand()%(SoluMulX_num+CandiNum);
    
    a11=RandIndex1;
    for(i=0;i<XNum;i++){
        if(a11<SoluMulX_new[i][0])
            break;
        a11-=SoluMulX_new[i][0];
    }
    a10=i;
    
    b11=RandIndex2;
    if(b11<SoluMulX_num){
        for(i=0;i<XNum;i++){
            if(b11<SoluMulX_new[i][0])
                break;
            b11-=SoluMulX_new[i][0];
        }
        b10=i;
        SwapInt(&SoluMulX_new[a10][1+a11],&SoluMulX_new[b10][1+b11]);
        r[0]=a10;
        r[1]=1+a11;
        r[2]=b10;
        r[3]=1+b11;
    }else{
        b11-=SoluMulX_num;
        SwapInt(&SoluMulX_new[a10][1+a11],&CandiSeq[b11]);
        r[0]=a10;
        r[1]=1+a11;
        r[2]=-1;
        r[3]=b11;
    }
}

void SA_MulXAccept(int** SoluMulX_cur, int** SoluMulX_new, const int XNum, double* E0, const double* E1){
    int i;
    for(i=0;i<XNum;i++){
        CopyInt(SoluMulX_new[i],SoluMulX_cur[i],SoluMulX_new[i][0]+1);
    }
    *E0=*E1;
}

void SA_aStep_MulX_NoCandiSeq(int** SoluMulX_cur, const int XNum, const double* MulX_DInf, double** SpotsPos, const int SpotsNum, double** MulX_PosCab, double* PTD_E0_Add, const double Temperature, int** SoluMulX_new, double* DistVec, double* RouteProb){
    int i,AcceptFlag,rr[2];
    double PTD_E1,DeltaE;
    
    for(i=0;i<XNum;i++){
        CopyInt(SoluMulX_cur[i],SoluMulX_new[i],SoluMulX_cur[i][0]+1);
    }
    SA_GloMove_MulXSwap_NoCandiSeq(SoluMulX_new,XNum,SpotsNum,rr);
    
    DeltaE=PTD_OneX(SoluMulX_new[rr[0]],MulX_DInf[rr[0]],SpotsPos,MulX_PosCab[rr[0]],DistVec,RouteProb)-PTD_OneX(SoluMulX_cur[rr[0]],MulX_DInf[rr[0]],SpotsPos,MulX_PosCab[rr[0]],DistVec,RouteProb);
    if(rr[1]>=0&&rr[1]!=rr[0]){
        DeltaE+=PTD_OneX(SoluMulX_new[rr[1]],MulX_DInf[rr[1]],SpotsPos,MulX_PosCab[rr[1]],DistVec,RouteProb)-PTD_OneX(SoluMulX_cur[rr[1]],MulX_DInf[rr[1]],SpotsPos,MulX_PosCab[rr[1]],DistVec,RouteProb);
    }
    
    PTD_E1=*PTD_E0_Add+DeltaE;
    AcceptFlag=SA_AcceptOrNot(*PTD_E0_Add,PTD_E1,Temperature);
    if(AcceptFlag==1){
        SA_MulXAccept(SoluMulX_cur,SoluMulX_new,XNum,PTD_E0_Add,&PTD_E1);//SoluMulX_cur, and PTD_E0 are the output of SA_aStep_MulX_NoCandiSeq.
        AccNum++;
    }
}

void SA_aStep_MulX(int** SoluMulX_cur, const int XNum, int* CandiSeq, const int CandiNum, const double* MulX_DInf, double** SpotsPos, double** MulX_PosCab, double* PTD_E0_Add, const double Temperature, int** SoluMulX_new, double* DistVec, double* RouteProb){
    int i,AcceptFlag,r[4];
    double PTD_E1,DeltaE;
    
    for(i=0;i<XNum;i++){
        CopyInt(SoluMulX_cur[i],SoluMulX_new[i],SoluMulX_cur[i][0]+1);
    }
    SA_GloMove_MulXSwap(SoluMulX_new,XNum,CandiSeq,CandiNum,r);
    
    DeltaE=PTD_OneX(SoluMulX_new[r[0]],MulX_DInf[r[0]],SpotsPos,MulX_PosCab[r[0]],DistVec,RouteProb)-PTD_OneX(SoluMulX_cur[r[0]],MulX_DInf[r[0]],SpotsPos,MulX_PosCab[r[0]],DistVec,RouteProb);
    if(r[2]>=0&&r[2]!=r[0]){
        DeltaE+=PTD_OneX(SoluMulX_new[r[2]],MulX_DInf[r[2]],SpotsPos,MulX_PosCab[r[2]],DistVec,RouteProb)-PTD_OneX(SoluMulX_cur[r[2]],MulX_DInf[r[2]],SpotsPos,MulX_PosCab[r[2]],DistVec,RouteProb);
    }
    
    PTD_E1=*PTD_E0_Add+DeltaE;
    AcceptFlag=SA_AcceptOrNot(*PTD_E0_Add,PTD_E1,Temperature);
    if(AcceptFlag==1){
        SA_MulXAccept(SoluMulX_cur,SoluMulX_new,XNum,PTD_E0_Add,&PTD_E1);//SoluMulX_cur, CandiSeq, and PTD_E0 are the output of SA_aStep_MulX.
        AccNum++;
    }else if(r[2]==-1){
        SwapInt(&SoluMulX_new[r[0]][r[1]],&CandiSeq[r[3]]);
    }
}

void SA_aStep_MulX_PMode2(int** SoluMulX_cur, const int XNum, int* CandiSeq, const int CandiNum, const double* MulX_DInf, double** SpotsPos, double** MulX_PosCab, double* PTD_E0_Add, const double Temperature, int** SoluMulX_new, double* DistVec, double* RouteProb, const int* OriSoluVec, int OriSoluLength, int** SoluMulXCandi_cur, int** SoluMulXCandi_new){
    int i,AcceptFlag,rr[2];
    double PTD_E1,DeltaE;
    
    for(i=0;i<XNum;i++){
        CopyInt(SoluMulX_cur[i],SoluMulX_new[i],SoluMulX_cur[i][0]+1);
        CopyInt(SoluMulXCandi_cur[i],SoluMulXCandi_new[i],SoluMulXCandi_cur[i][0]+1);
    }
    
    SA_GloMove_MulXSwap_PMode2(SoluMulX_new,XNum,CandiSeq,CandiNum,OriSoluVec,OriSoluLength,SoluMulXCandi_new,rr);
    
    DeltaE=PTD_OneX(SoluMulX_new[rr[0]],MulX_DInf[rr[0]],SpotsPos,MulX_PosCab[rr[0]],DistVec,RouteProb)-PTD_OneX(SoluMulX_cur[rr[0]],MulX_DInf[rr[0]],SpotsPos,MulX_PosCab[rr[0]],DistVec,RouteProb);
    if(rr[1]>=0&&rr[1]!=rr[0]){
        DeltaE+=PTD_OneX(SoluMulX_new[rr[1]],MulX_DInf[rr[1]],SpotsPos,MulX_PosCab[rr[1]],DistVec,RouteProb)-PTD_OneX(SoluMulX_cur[rr[1]],MulX_DInf[rr[1]],SpotsPos,MulX_PosCab[rr[1]],DistVec,RouteProb);
    }
    
    PTD_E1=*PTD_E0_Add+DeltaE;
    AcceptFlag=SA_AcceptOrNot(*PTD_E0_Add,PTD_E1,Temperature);
    if(AcceptFlag==1){
        SA_MulXAccept(SoluMulX_cur,SoluMulX_new,XNum,PTD_E0_Add,&PTD_E1);
        SA_MulXAccept(SoluMulXCandi_cur,SoluMulXCandi_new,XNum,PTD_E0_Add,&PTD_E1);//SoluMulX_cur, SoluMulXCandi_cur, CandiSeq, and PTD_E0 are the output of SA_aStep_MulX.
        AccNum++;
    }
}

// Function Series 4: Parallel
// Function Series 4.1: My Techniques

void GRoute_GDInf_GCabsPos_update(int** GRoute, double** GCabsPos, double* GDInf, int** Route, double** CabsPos, const double* DInf, const int CurRouteIndex, const int CabsNumPerCore){
    int i,j;
    
    j=0;
    for(i=CurRouteIndex;i<CurRouteIndex+CabsNumPerCore;i++){
        GRoute[j][0]=Route[i][0];
        GCabsPos[j][0]=CabsPos[i][0];
        GCabsPos[j][1]=CabsPos[i][1];
        GDInf[j]=DInf[i];
        j++;
    }
}

void Parallel_rotation(int* Route, int* TempRoute, const int GRank, const int GSize, const int Mode_ISend, MPI_Comm GroupWorld){
    
    int i,RouteVecLength;
    MPI_Request send_request,recv_request;
    RouteVecLength=RouteLimit*CabsNum/GSize;
    
    if(GSize<2)
        printf("Error! GSize should be >=2 to use the function Parallel_rotation!\n");
    
    if(Mode_ISend==0){
        if (GRank==0){
            MPI_Recv(TempRoute,RouteVecLength,MPI_INT,GSize-1,MPI_ANY_TAG,GroupWorld,MPI_STATUS_IGNORE);
            MPI_Send(Route,RouteVecLength,MPI_INT,GRank+1,0,GroupWorld);
        }
        else if(GRank==GSize-1){
            MPI_Send(Route,RouteVecLength,MPI_INT,0,0,GroupWorld);
            MPI_Recv(TempRoute,RouteVecLength,MPI_INT,GRank-1,MPI_ANY_TAG,GroupWorld,MPI_STATUS_IGNORE);
        }
        else{
            MPI_Recv(TempRoute,RouteVecLength,MPI_INT,GRank-1,MPI_ANY_TAG,GroupWorld,MPI_STATUS_IGNORE);
            MPI_Send(Route,RouteVecLength,MPI_INT,GRank+1,0,GroupWorld);
        }
    }else if(Mode_ISend==1){
        if(GRank==0){
            MPI_Isend(Route,RouteVecLength,MPI_INT,GRank+1,0,GroupWorld,&send_request);
            MPI_Irecv(TempRoute,RouteVecLength,MPI_INT,GSize-1,MPI_ANY_TAG,GroupWorld,&recv_request);
        }
        else if(GRank==GSize-1){
            MPI_Isend(Route,RouteVecLength,MPI_INT,0,0,GroupWorld,&send_request);
            MPI_Irecv(TempRoute,RouteVecLength,MPI_INT,GRank-1,MPI_ANY_TAG,GroupWorld,&recv_request);
        }
        else{
            MPI_Isend(Route,RouteVecLength,MPI_INT,GRank+1,0,GroupWorld,&send_request);
            MPI_Irecv(TempRoute,RouteVecLength,MPI_INT,GRank-1,MPI_ANY_TAG,GroupWorld,&recv_request);
        }
        MPI_Wait(&send_request,MPI_STATUS_IGNORE);
        MPI_Wait(&recv_request,MPI_STATUS_IGNORE);
    }else{
        printf("Error: Mode_ISend is neither 0 or 1.\n");
    }
    for(i=0;i<RouteVecLength;i++)
        Route[i]=TempRoute[i];
}

void KnuthShuffles(int* Seq, const int n){
    int i,temp;
    for(i=0;i<n-1;i++){
        temp=rand()%(n-i);
        SwapInt(&Seq[i],&Seq[i+temp]);
    }
}

// Route length should be the same! RouteLengthPerCab=Route[i][0]+1.
void Parallel_Interact(int** GRoute, const int CabsNumPerCore, const int RouteLengthPerCab, int ShiftTimes, int* RouteVec, int* RouteVecTemp, int* TempVec, int* TempVec2, double** CabsPos, int** Route, int** RouteNew, const double* DInf, double** SpotsPos, double* PTD_E0_Add, const int StepsPerInter, const int RotateOrNot, const int GRank, const int GSize, MPI_Comm GroupWorld, double* DistVec, double* RouteProb, double* TemperatureAdd, double PCoolRate, const int Mode_ParaTech, const int InGRank, const int LGRank, const int InGroupSize, MPI_Comm InGWorld, MPI_Comm LGWorld, int* OptRankAdd){
    int step,i,j,k,Length,*CandiSeq,Index;
    struct pair_ERank PIn;
    struct pair_ERank POut;
    
    ShiftTimes-=GSize;
    if(ShiftTimes>=0){
        printf("Error! In Parallel_Interact, ShiftTimes>=0!\n");
        fflush(stdout);
    }
    
    Length=RouteLengthPerCab*CabsNumPerCore;
    Route2RouteVec(GRoute,RouteVec,CabsNumPerCore);
    MPI_Allgather(RouteVec,Length,MPI_INT,RouteVecTemp,Length,MPI_INT,GroupWorld);
    for(i=0;i<GSize;i++){
        Index=((i-ShiftTimes)%GSize)*CabsNumPerCore*RouteLengthPerCab;
        for(j=0;j<CabsNumPerCore;j++){
            for(k=0;k<RouteLengthPerCab;k++){
                Route[i*CabsNumPerCore+j][k]=RouteVecTemp[Index+j*RouteLengthPerCab+k];
            }
        }
    }
    
    *PTD_E0_Add=PTD_MulX(Route,DInf,SpotsPos,CabsPos,CabsNum,DistVec,RouteProb);
    for(step=0;step<StepsPerInter;step++){
        SA_aStep_MulX(Route,CabsNum,CandiSeq,0,DInf,SpotsPos,CabsPos,PTD_E0_Add,*TemperatureAdd,RouteNew,DistVec,RouteProb);
        *TemperatureAdd=(*TemperatureAdd)*PCoolRate;
    }
    
    if(RotateOrNot==1)
        ShiftTimes--;
    for(i=0;i<GSize;i++){
        Index=((i-ShiftTimes)%GSize)*CabsNumPerCore*RouteLengthPerCab;
        for(j=0;j<CabsNumPerCore;j++){
            for(k=0;k<RouteLengthPerCab;k++){
                RouteVecTemp[Index+j*RouteLengthPerCab+k]=Route[i*CabsNumPerCore+j][k];
            }
        }
    }
    
    PIn.E=*PTD_E0_Add;
    
    if(Mode_ParaTech==1){
        PIn.rank=GRank;
        MPI_Allreduce(&PIn,&POut,1,MPI_DOUBLE_INT,MPI_MINLOC,GroupWorld);
    }else if(Mode_ParaTech==2){
        PIn.rank=LGRank;
        MPI_Allreduce(&PIn,&POut,1,MPI_DOUBLE_INT,MPI_MINLOC,LGWorld);
    }
    
    if(Mode_ParaTech==2&&InGRank==POut.rank%InGroupSize&&GRank==POut.rank/InGroupSize){
        if(LGRank!=POut.rank)
            printf("Error! In Parallel_Interact, LGRank, GRank, and InGRank are wrong!\n");
        if(*PTD_E0_Add!=POut.E)
            printf("Error! In Parallel_Interact, the optimal PTD_E0 should == POut.E but it is not!\n");
    }
    *PTD_E0_Add=POut.E;
    
    if(Mode_ParaTech==1){
        MPI_Scatter(RouteVecTemp,Length,MPI_INT,RouteVec,Length,MPI_INT,POut.rank,GroupWorld);
    }else if(Mode_ParaTech==2){
        if(InGRank==POut.rank%InGroupSize)
            MPI_Scatter(RouteVecTemp,Length,MPI_INT,RouteVec,Length,MPI_INT,POut.rank/InGroupSize,GroupWorld);
        MPI_Bcast(RouteVec,Length,MPI_INT,POut.rank%InGroupSize,InGWorld);
    }
    
    RouteVec2Route(GRoute,RouteVec,CabsNumPerCore);
    *OptRankAdd=POut.rank;
}

// Load imbalance problem is in this function.
// GRank0_StoreCandiSeq malloc length should be total number of pick-up points. TempVec and TempVec2 have the size at least GSize.
void Parallel_Shuffle(int* CandiSeq, int CandiNum, int* GRank0_StoreCandiSeq, int* TempVec, int* TempVec2, const int GRank, const int GSize, MPI_Comm GroupWorld, const int Mode_ParaTech, const int InGRank, MPI_Comm InGWorld){
    
    int i,GRank0_StoreCandiSeqLength;
    
    if(Mode_ParaTech==1||Mode_ParaTech==2&&InGRank==0){
        
        MPI_Gather(&CandiNum,1,MPI_INT,TempVec,1,MPI_INT,0,GroupWorld);
        if(GRank==0){
            TempVec2[0]=0;
            for(i=0;i<GSize-1;i++){
                TempVec2[i+1]=TempVec2[i]+TempVec[i];
            }
        }
        MPI_Gatherv(CandiSeq,CandiNum,MPI_INT,GRank0_StoreCandiSeq,TempVec,TempVec2,MPI_INT,0,GroupWorld);
        if(GRank==0){
            GRank0_StoreCandiSeqLength=TempVec[GSize-1]+TempVec2[GSize-1];
            
            KnuthShuffles(GRank0_StoreCandiSeq,GRank0_StoreCandiSeqLength);
            
            CandiNum=GRank0_StoreCandiSeqLength/GSize;
            TempVec2[0]=0;
            for(i=0;i<GRank0_StoreCandiSeqLength%GSize;i++){
                TempVec[i]=CandiNum+1;
                TempVec2[i+1]=TempVec2[i]+TempVec[i];// Note that i will not be GSize-1.
            }
            for(i=GRank0_StoreCandiSeqLength%GSize;i<GSize-1;i++){
                TempVec[i]=CandiNum;
                TempVec2[i+1]=TempVec2[i]+TempVec[i];// Note that i will not be GSize-1.
            }
            TempVec[GSize-1]=CandiNum;
        }
        
        MPI_Scatter(TempVec,1,MPI_INT,&CandiNum,1,MPI_INT,0,GroupWorld);
        MPI_Scatterv(GRank0_StoreCandiSeq,TempVec,TempVec2,MPI_INT,CandiSeq,CandiNum,MPI_INT,0,GroupWorld);
        
    }
    
    if(Mode_ParaTech==2){
        
        MPI_Bcast(CandiSeq,CandiNum,MPI_INT,0,InGWorld);
        
    }
}

// Function Series 4.2: My Mixing Period and Mixing Pattern

void GRouteCandi_OriSoluVec_Ini(int** GRouteCandi, int* OriSoluVec, int** GRoute, int CabsNumPerCore){
    int i,j,k;
    
    for(i=0,k=0;i<CabsNumPerCore;i++){
        GRouteCandi[i][0]=GRoute[i][0];
        for(j=1;j<=GRoute[i][0];j++){
            OriSoluVec[k]=GRoute[i][j];
            GRouteCandi[i][j]=k++;
        }
    }
}

// SearchSteps is a parameter.
double Parallel_SimilarityRate_OneX(const int X_dim, const int CandiNum, const double AccRate, const int SearchSteps){
    double temp,l,n,r,m;
    
    l=(double)X_dim;
    n=(double)CandiNum;
    r=1.0-AccRate;
    m=(double)SearchSteps;
    
    temp=(l-1)/l*((n-1)/n+1/n*r)+(1/l)*r;
    return pow(temp,m);
}

// X_temp is at least X_dim.
void Parallel_SimilarityMPattern_OneX(int* X, int* X_temp, const int X_dim, double* EAdd, double SRate, int* ReplaceOrNotAdd, const int GRank, MPI_Comm MPIWorld){
    double TwoVecSRate;
    struct pair_ERank PIn;
    struct pair_ERank POut;
    
    CopyInt(X,X_temp,X_dim);
    PIn.E=*EAdd;
    PIn.rank=GRank;
    MPI_Allreduce(&PIn,&POut,1,MPI_DOUBLE_INT,MPI_MINLOC,MPIWorld);
    MPI_Bcast(&SRate,1,MPI_DOUBLE,POut.rank,MPIWorld);
    MPI_Bcast(X_temp,X_dim,MPI_INT,POut.rank,MPIWorld);
    TwoVecSRate=TwoVec_SimilarityRate(X,X_temp,X_dim);
    if(TwoVecSRate<SRate){
        CopyInt(X_temp,X,X_dim);
        *EAdd=POut.E;
        *ReplaceOrNotAdd=1;
    }else{
        *ReplaceOrNotAdd=0;
    }
    if(POut.rank==GRank&&*ReplaceOrNotAdd!=0){
        printf("Error! In Parallel_SimilarityMPattern_OneX, my understanding of ReplaceOrNot is wrong.\n");
        fflush(stdout);
    }
}

// ReplaceInfoVec is at least GSize. InRangeRate is a parameter.
void Parallel_SimilarityMPeriod_OneX(const int ReplaceOrNot, int* ReplaceInfoVec, const double InRangeRate, int* MPeriodAdd, const int GSize, MPI_Comm MPIWorld){
    int i,sum;
    
    MPI_Allgather(&ReplaceOrNot,1,MPI_INT,ReplaceInfoVec,1,MPI_INT,MPIWorld);
    for(i=0,sum=-1;i<GSize;i++){
        sum+=1-ReplaceInfoVec[i];
    }
    if(((double)sum)/(double)(GSize-1)<InRangeRate)
        (*MPeriodAdd)--;
    else
        (*MPeriodAdd)++;
    if(*MPeriodAdd<1)
        *MPeriodAdd=1;
}

void Parallel_OneToAll(int* X, const int X_dim, double* EAdd, const int GRank, MPI_Comm MPIWorld){
    struct pair_ERank PIn;
    struct pair_ERank POut;
    
    PIn.E=*EAdd;
    PIn.rank=GRank;
    MPI_Allreduce(&PIn,&POut,1,MPI_DOUBLE_INT,MPI_MINLOC,MPIWorld);
    MPI_Bcast(X,X_dim, MPI_INT,POut.rank,MPIWorld);
    *EAdd=POut.E;
}

// Function Series 4.3: CDR
double CDR_ComputeP(const double *a, const int l, const int r, const double t){
    int i;
    double p,m;
    p=0.0;
    if (t>0.00000001){
        for(i=0;i<l;i++)
            p+=exp((a[r]-a[i])/t);
        return 1.0/p;
    }else{
        m=a[0];
        for(i=1;i<l;i++)
            if(a[i]<m)
                m=a[i];
        if(a[r]==m){
            for(i=0;i<l;i++)
                if(a[i]==m)
                    p+=1.0;
            return 1.0/p;
        }else
            return 0.0;
    }
}

void Parallel_CDR(int* X, int* X_AllCores, const int X_dim, double* EAdd, double* E_AllCores, double* EProb_AllCores, const double temperature, int* ChosenStateAdd, const int GRank, const int GSize, MPI_Comm MPIWorld){
    int i,j,temp2;
    double EProb,EProb_sum,temp1;
    
    MPI_Allgather(X,X_dim,MPI_INT,X_AllCores,X_dim,MPI_INT,MPIWorld);
    MPI_Allgather(EAdd,1,MPI_DOUBLE,E_AllCores,1,MPI_DOUBLE,MPIWorld);
    
    EProb=CDR_ComputeP(E_AllCores,GSize,GRank,temperature);
    MPI_Allgather(&EProb,1,MPI_DOUBLE,EProb_AllCores,1,MPI_DOUBLE,MPIWorld);
    temp1=((double)rand())/RAND_MAX;
    EProb_sum=0.0;
    for(i=0;i<GSize;i++){
        EProb_sum+=EProb_AllCores[i];
        if(temp1<=EProb_sum)
            break;
    }
    if(i>=GSize){
        printf("Error! CDR communication is wrong.\n");
        fflush(stdout);
    }
    *EAdd=E_AllCores[i];
    temp2=i*X_dim;
    for(j=temp2;j<(i+1)*X_dim;j++)
        X[j-temp2]=X_AllCores[j];
    *ChosenStateAdd=i;
}

// Function Series 4.4: Lou
// StatesVec should be at least GSize.
void Parallel_RI(int* MPeriodAdd, const int ChosenState, int* StatesVec, const int MaxMP, const int GSize, MPI_Comm MPIWorld){
    int i,distinct,CurVal;
    double Rnext,Rcur,C,Ar,A0;
    //log(Rnext)=log(Rcur)+C(Ar-A0);
    
    C=2*log(2);
    A0=0.5;
    
    MPI_Allgather(&ChosenState,1,MPI_INT,StatesVec,1,MPI_INT,MPIWorld);
    SortInt_NoPlace(StatesVec,GSize);
    distinct=1;
    CurVal=StatesVec[0];
    for(i=1;i<GSize;i++){
        if(CurVal!=StatesVec[i]){
            distinct++;
            CurVal=StatesVec[i];
        }
    }
    Ar=((double)distinct)/(double)GSize;
    Rcur=(double)(*MPeriodAdd);
    Rnext=exp(log(Rcur)+C*(Ar-A0));
    *MPeriodAdd=(int)round(Rnext);
    if(*MPeriodAdd<1)
        *MPeriodAdd=1;
    else if(*MPeriodAdd>MaxMP)
        *MPeriodAdd=MaxMP;
}

// Function Series 5: Seq_only, CDR_Lou, and MyTech_MyMPs.
void Seq_only(double** CabsPos, int** Route, const double* DInf, double** SpotsPos, const int SpotsNum, const int Mode_PrintProc, const int PrintPerSteps, const int StopNum, double Temperature, const double Alpha, FILE* ProcFile, double* PTD_E0_Add, double* PTD_min_Add, int* StepNum_Add){
    int i,**RouteNew;
    double PTD_Record,DistVec[RouteLimit],RouteProb[RouteLimit];
    
    int StepNum,CandiNum,*CandiSeq,StopStep;
    double PTD_E0,PTD_min;
    
    RouteNew=(int**)malloc(CabsNum*sizeof(int*));
    for(i=0;i<CabsNum;i++){
        RouteNew[i]=(int*)malloc(RouteLimit*sizeof(int));
    }
    CandiSeq=(int*)malloc(SpotsNum*sizeof(int));
    
    for(i=0;i<SpotsNum;i++)
        CandiSeq[i]=i;
    KnuthShuffles(CandiSeq,SpotsNum);
    CandiNum=SpotsNum;
    
    SA_MulXIni(Route,CabsNum,CandiSeq,&CandiNum);
    PTD_E0=PTD_MulX(Route,DInf,SpotsPos,CabsPos,CabsNum,DistVec,RouteProb);
    PTD_min=1000000.0;
    StepNum=0;
    AccNum=0;
    StopStep=0;
    PTD_Record=0.0;
    
    while(1){
        StepNum++;
        SA_aStep_MulX(Route,CabsNum,CandiSeq,CandiNum,DInf,SpotsPos,CabsPos,&PTD_E0,Temperature,RouteNew,DistVec,RouteProb);
        
        if(Mode_PrintProc==1){
            if(StepNum%PrintPerSteps==1){
                fprintf(ProcFile,"%f %d\n",PTD_E0,StepNum);
                fflush(ProcFile);
            }
        }
        
        Temperature=Temperature*Alpha;
        
        if(PTD_min>PTD_E0)
            PTD_min=PTD_E0;
        
        if(PTD_E0==PTD_Record)
            StopStep++;
        else{
            StopStep=0;
            PTD_Record=PTD_E0;
        }
        if(StopStep>=StopNum)
            break;
    }
    
    *PTD_E0_Add=PTD_E0;
    *PTD_min_Add=PTD_min;
    *StepNum_Add=StepNum;
    
    free(CandiSeq);
    for(i=0;i<CabsNum;i++){
        free(RouteNew[i]);
    }
    free(RouteNew);
}

void CDR_Lou(double** CabsPos, int** Route, const double* DInf, double** SpotsPos, const int SpotsNum, const int Mode_PrintProc, const int PrintPerSteps, const int StopNum, double Temperature, const double Alpha, FILE* ProcFile, double* PTD_E0_Add, double* PTD_min_Add, int* StepNum_Add, const int Mode_ParaTech, int MPeriod, const int MaxMP, const int GRank, const int GroupSize, MPI_Comm GWorld){
    int i,j,ChosenState,temp1,tempI,**RouteNew;
    double PTD_Record,tempE,DistVec[RouteLimit],RouteProb[RouteLimit];
    
    int StepNum,StopStep,*RouteVec,RouteVecLength,*RouteVecAllCores,*StatesVec;
    double PTD_E0,PTD_min,PCoolRate,*E_AllCores,*EProb_AllCores;
    
    for(i=0,RouteVecLength=0;i<CabsNum;i++){
        RouteVecLength+=1+Route[i][0];
    }
    
    RouteNew=(int**)malloc(CabsNum*sizeof(int*));
    for(i=0;i<CabsNum;i++){
        RouteNew[i]=(int*)malloc(RouteLimit*sizeof(int));
    }
    RouteVec=(int*)malloc(RouteVecLength*sizeof(int));
    RouteVecAllCores=(int*)malloc(GroupSize*RouteVecLength*sizeof(int));
    E_AllCores=(double*)malloc(GroupSize*sizeof(double));
    EProb_AllCores=(double*)malloc(GroupSize*sizeof(double));
    StatesVec=(int*)malloc(GroupSize*sizeof(int));
    
    SA_MulXIni_NoCandiSeq(Route,CabsNum,SpotsNum);
    PTD_E0=PTD_MulX(Route,DInf,SpotsPos,CabsPos,CabsNum,DistVec,RouteProb);
    PTD_min=1000000.0;
    StepNum=0;
    AccNum=0;
    StopStep=0;
    PTD_Record=0.0;
    PCoolRate=pow(Alpha,GroupSize);
    
    while(1){
        StepNum++;
        SA_aStep_MulX_NoCandiSeq(Route,CabsNum,DInf,SpotsPos,SpotsNum,CabsPos,&PTD_E0,Temperature,RouteNew,DistVec,RouteProb);
        
        if(GRank==0&&Mode_PrintProc==1){
            if(StepNum%PrintPerSteps==1){
                fprintf(ProcFile,"%f %d\n",PTD_E0,StepNum);
                fflush(ProcFile);
            }
        }
        
        Temperature=Temperature*PCoolRate;
        
        if(StepNum%MPeriod==0){
            
            Route2RouteVec(Route,RouteVec,CabsNum);
            Parallel_CDR(RouteVec,RouteVecAllCores,RouteVecLength,&PTD_E0,E_AllCores,EProb_AllCores,Temperature,&ChosenState,GRank,GroupSize,GWorld);
            RouteVec2Route(Route,RouteVec,CabsNum);
            
            if(Mode_ParaTech==4)
                Parallel_RI(&MPeriod,ChosenState,StatesVec,MaxMP,GroupSize,GWorld);
            
            tempE=E_AllCores[0];
            tempI=0;
            for(i=1;i<GroupSize;i++){
                if(E_AllCores[i]<tempE){
                    tempE=E_AllCores[i];
                    tempI=i;
                }
                    
            }
            if(PTD_min>tempE)
                PTD_min=tempE;
            
            if(tempE==PTD_Record)
                StopStep++;
            else{
                StopStep=0;
                PTD_Record=tempE;
            }
            if(StopStep*GroupSize*MPeriod>=StopNum)
                break;
        }
    }
    
    temp1=tempI*RouteVecLength;
    for(j=temp1;j<(tempI+1)*RouteVecLength;j++)
        RouteVec[j-temp1]=RouteVecAllCores[j];
    RouteVec2Route(Route,RouteVec,CabsNum);
    
    *PTD_E0_Add=tempE;
    *PTD_min_Add=PTD_min;
    *StepNum_Add=StepNum;
    
    free(RouteVec);
    free(RouteVecAllCores);
    free(E_AllCores);
    free(EProb_AllCores);
    free(StatesVec);
    for(i=0;i<CabsNum;i++){
        free(RouteNew[i]);
    }
    free(RouteNew);
}

// AccRate here is initial AccRate.
// Although each core only computes part of the routes, they claim the memory for all route.
void MyTech_MyMPs(double** CabsPos, int** Route, const double* DInf, double** SpotsPos, const int SpotsNum, const int Mode_PrintProc, const int PrintPerSteps, const int StopNum, double Temperature, const double Alpha, FILE* ProcFile, double* PTD_E0_Add, double* PTD_min_Add, int* StepNum_Add, const int Mode_ParaTech, const int GRank, const int GroupSize, MPI_Comm GWorld, const int LGRank, const int LGroupSize, MPI_Comm LGWorld, const int InGRank, const int InGroupSize, MPI_Comm InGWorld, const int Mode_ISend, const int Steps_PerRo, const int RoTimes_PerInter, const int RoTimes_PerShulffle, const int SASteps_InInter, const int SearchSteps, const double InRangeRate, int MPeriod, const int AccPeriod, double AccRate, int* OptRankAdd){
    
    int i,j,k,*TempVec,*TempVec2,**GRouteNew,*RouteVec,*RouteVecTemp,RouteVecLength,**RouteNew,FlagTemp;
    double PTD_Record,DistVec[RouteLimit],RouteProb[RouteLimit],*GDInf,**GCabsPos;
    
    int StepNum,CandiNum,*GRank0_StoreCandiSeq,*CandiSeq,StopStep,CabsNumPerCore,**GRoute,CurRouteIndex,ShiftTimes,NeedRoOrNOt,ShiftOrNot,RoTimes;
    double PTD_E0,PTD_min,PCoolRate,PTD_E0_AllRoute;
    
    int *OriSoluVec,OriSoluLength,**GRouteCandi,**GRouteCandi_new,ReplaceOrNot,*ReplaceInfoVec,*XVec1,*XVec2,*XVec3,Len;
    double SRate;
    
    CabsNumPerCore=CabsNum/GroupSize;
    ShiftTimes=0;
    CurRouteIndex=CabsNumPerCore*((GRank+ShiftTimes)%GroupSize);
    for(i=0,RouteVecLength=0;i<CabsNum;i++)
        RouteVecLength+=1+Route[i][0];
    if(Mode_ParaTech==2)
        OriSoluLength=CabsNumPerCore*Route[0][0];// All cabs route lengths are the same.
    
    if(GRank==0){
        GRank0_StoreCandiSeq=(int*)malloc(SpotsNum*sizeof(int));
    }
    RouteNew=(int**)malloc(CabsNum*sizeof(int*));
    for(i=0;i<CabsNum;i++){
        RouteNew[i]=(int*)malloc(RouteLimit*sizeof(int));
    }
    CandiSeq=(int*)malloc((SpotsNum/GroupSize+1)*sizeof(int));
    TempVec=(int*)malloc(GroupSize*sizeof(int));
    TempVec2=(int*)malloc(GroupSize*sizeof(int));
    GRouteNew=(int**)malloc(CabsNumPerCore*sizeof(int*));
    for(i=0;i<CabsNumPerCore;i++){
        GRouteNew[i]=(int*)malloc(RouteLimit*sizeof(int));
    }
    GRoute=(int**)malloc(CabsNumPerCore*sizeof(int*));
    for(i=0;i<CabsNumPerCore;i++){
        GRoute[i]=(int*)malloc(RouteLimit*sizeof(int));
    }
    GDInf=(double*)malloc(CabsNumPerCore*sizeof(double));
    GCabsPos=(double**)malloc(CabsNumPerCore*sizeof(double*));
    for(i=0;i<CabsNumPerCore;i++){
        GCabsPos[i]=(double*)malloc(2*sizeof(double));
    }
    RouteVec=(int*)malloc(RouteVecLength*sizeof(int));
    RouteVecTemp=(int*)malloc(RouteVecLength*sizeof(int));
    OriSoluVec=(int*)malloc(OriSoluLength*sizeof(int));
    XVec1=(int*)malloc(OriSoluLength*2*sizeof(int));
    XVec2=(int*)malloc(OriSoluLength*2*sizeof(int));
    ReplaceInfoVec=(int*)malloc(InGroupSize*sizeof(int));
    GRouteCandi=(int**)malloc(CabsNumPerCore*sizeof(int*));
    for(i=0;i<CabsNumPerCore;i++){
        GRouteCandi[i]=(int*)malloc(RouteLimit*sizeof(int));
    }
    GRouteCandi_new=(int**)malloc(CabsNumPerCore*sizeof(int*));
    for(i=0;i<CabsNumPerCore;i++){
        GRouteCandi_new[i]=(int*)malloc(RouteLimit*sizeof(int));
    }
    
    FlagTemp=0;
    if(Mode_ParaTech==1&&GRank==0)
        FlagTemp=1;
    else if(Mode_ParaTech==2&&LGRank==0)
        FlagTemp=1;
    
    if(FlagTemp==1){
        for(i=0;i<SpotsNum;i++)
            GRank0_StoreCandiSeq[i]=i;
        KnuthShuffles(GRank0_StoreCandiSeq,SpotsNum);
        
        CandiNum=SpotsNum/GroupSize;
        TempVec2[0]=0;
        for(i=0;i<SpotsNum%GroupSize;i++){
            TempVec[i]=CandiNum+1;
            TempVec2[i+1]=TempVec2[i]+TempVec[i];// Note that i will not be GroupSize-1.
        }
        for(i=SpotsNum%GroupSize;i<GroupSize-1;i++){
            TempVec[i]=CandiNum;
            TempVec2[i+1]=TempVec2[i]+TempVec[i];// Note that i will not be GroupSize-1.
        }
        TempVec[GroupSize-1]=CandiNum;
    }
    CandiNum=SpotsNum/GroupSize;
    if(GRank<SpotsNum%GroupSize)
        CandiNum++;
    
    FlagTemp=0;
    if(Mode_ParaTech==1)
        FlagTemp=1;
    else if(Mode_ParaTech==2&&InGRank==0)
        FlagTemp=1;
    
    if(FlagTemp==1)
        MPI_Scatterv(GRank0_StoreCandiSeq,TempVec,TempVec2,MPI_INT,CandiSeq,CandiNum,MPI_INT,0,GWorld);
    if(Mode_ParaTech==2)
        MPI_Bcast(CandiSeq,CandiNum,MPI_INT,0,InGWorld);
    
    GRoute_GDInf_GCabsPos_update(GRoute,GCabsPos,GDInf,Route,CabsPos,DInf,CurRouteIndex,CabsNumPerCore);

    for(i=0;i<CabsNumPerCore;i++){
        for(j=1;j<=GRoute[i][0];j++)
            GRoute[i][j]=CandiSeq[-1+CandiNum--];
    }
    
    if(Mode_ParaTech==2)
        GRouteCandi_OriSoluVec_Ini(GRouteCandi,OriSoluVec,GRoute,CabsNumPerCore);
    
    PTD_E0=PTD_MulX(GRoute,GDInf,SpotsPos,GCabsPos,CabsNumPerCore,DistVec,RouteProb);
    PTD_min=1000000.0;
    StepNum=0;
    AccNum=0;
    StopStep=0;
    PTD_Record=0.0;
    if(Mode_ParaTech==1)
        PCoolRate=pow(Alpha,GroupSize);
    else if(Mode_ParaTech==2)
        PCoolRate=pow(Alpha,LGroupSize);
    RoTimes=0;
    ShiftTimes=0;// 0 is the only common one that all GRanks can reach.
    
    while(1){
        StepNum++;
        if(Mode_ParaTech==1)
            SA_aStep_MulX(GRoute,CabsNumPerCore,CandiSeq,CandiNum,GDInf,SpotsPos,GCabsPos,&PTD_E0,Temperature,GRouteNew,DistVec,RouteProb);
        else if(Mode_ParaTech==2)
            SA_aStep_MulX_PMode2(GRoute,CabsNumPerCore,CandiSeq,CandiNum,GDInf,SpotsPos,GCabsPos,&PTD_E0,Temperature,GRouteNew,DistVec,RouteProb,OriSoluVec,OriSoluLength,GRouteCandi,GRouteCandi_new);
        
        if(GRank==0&&Mode_PrintProc==1){
            if(StepNum%PrintPerSteps==1){
                if(Mode_ParaTech==1){
                    fprintf(ProcFile,"%f %d\n",PTD_E0,StepNum);
                    fflush(ProcFile);
                }else if(Mode_ParaTech==2){
                    fprintf(ProcFile,"%f %d %f %d\n",PTD_E0,StepNum,AccRate,MPeriod);
                    fflush(ProcFile);
                }
            }
        }
        
        Temperature=Temperature*PCoolRate;
        
        if(Mode_ParaTech==2&&StepNum%AccPeriod==0){
            AccRate=((double)AccNum)/(double)AccPeriod;
            AccNum=0;
        }
        
        if(Mode_ParaTech==2&&StepNum%MPeriod==0){
            
            // ??? SRate suitable for multiple Xs. ??? Initial AccRate.
            SRate=Parallel_SimilarityRate_OneX(OriSoluLength,CandiNum,AccRate,SearchSteps);
            
            for(i=0;i<CabsNumPerCore;i++)
                CopyInt(&GRoute[i][1],&XVec1[i*GRoute[0][0]],GRoute[0][0]);
            for(i=0;i<CabsNumPerCore;i++)
                CopyInt(&GRouteCandi[i][1],&XVec1[OriSoluLength+i*GRoute[0][0]],GRoute[0][0]);
            
            Parallel_SimilarityMPattern_OneX(XVec1,XVec2,OriSoluLength*2,&PTD_E0,SRate,&ReplaceOrNot,InGRank,InGWorld);
            
            for(i=0;i<CabsNumPerCore;i++)
                CopyInt(&XVec1[i*GRoute[0][0]],&GRoute[i][1],GRoute[0][0]);
            for(i=0;i<CabsNumPerCore;i++)
                CopyInt(&XVec1[OriSoluLength+i*GRoute[0][0]],&GRouteCandi[i][1],GRoute[0][0]);
            
            Parallel_SimilarityMPeriod_OneX(ReplaceOrNot,ReplaceInfoVec,InRangeRate,&MPeriod,InGroupSize,InGWorld);
        }
        
        if(StepNum%Steps_PerRo==0){
            ShiftOrNot=0;
            NeedRoOrNOt=1;
            RoTimes++;
            
            if(Mode_ParaTech==2){
                
                for(i=0;i<CabsNumPerCore;i++)
                    CopyInt(&GRoute[i][1],&XVec1[i*GRoute[0][0]],GRoute[0][0]);
                for(i=0;i<CabsNumPerCore;i++)
                    CopyInt(&GRouteCandi[i][1],&XVec1[OriSoluLength+i*GRoute[0][0]],GRoute[0][0]);
                
                Parallel_OneToAll(XVec1,OriSoluLength*2,&PTD_E0,InGRank,InGWorld);
                
                for(i=0;i<CabsNumPerCore;i++)
                    CopyInt(&XVec1[i*GRoute[0][0]],&GRoute[i][1],GRoute[0][0]);
                
                XVec3=&XVec1[OriSoluLength];
                SortInt_NoPlace(XVec3,OriSoluLength);
                
                Len=OriSoluLength;
                for(i=0,j=0;i<OriSoluLength;i++){
                    if(i==XVec3[j]){
                        j++;
                    }else{
                        SwapInt(&OriSoluVec[i],&CandiSeq[XVec3[--Len]-OriSoluLength]);
                    }
                }
                if(j!=Len){
                    printf("Error! In Mode 2 main function, synchronize procedure is not as expected!\n");
                    fflush(stdout);
                }
            }
            
            if(RoTimes%RoTimes_PerInter==0){
                
                Parallel_Interact(GRoute,CabsNumPerCore,Route[0][0]+1,ShiftTimes,RouteVec,RouteVecTemp,TempVec,TempVec2,CabsPos,Route,RouteNew,DInf,SpotsPos,&PTD_E0_AllRoute,SASteps_InInter,1,GRank,GroupSize,GWorld,DistVec,RouteProb,&Temperature,PCoolRate,Mode_ParaTech,InGRank,LGRank,InGroupSize,InGWorld,LGWorld,OptRankAdd);
                
                //StepNum+=SASteps_InInter;// !!!!! Every step does not have the same amount of one cab PTD evaluation numbers.
                
                if(PTD_min>PTD_E0_AllRoute)
                    PTD_min=PTD_E0_AllRoute;
                
                if(PTD_E0_AllRoute==PTD_Record)
                    StopStep++;
                else{
                    StopStep=0;
                    PTD_Record=PTD_E0_AllRoute;
                }
                if(Mode_ParaTech==1){
                    if(StopStep*GroupSize*(Steps_PerRo*RoTimes_PerInter+SASteps_InInter)>=StopNum)
                        break;
                }else if(Mode_ParaTech==2){
                    if(StopStep*LGroupSize*(Steps_PerRo*RoTimes_PerInter+SASteps_InInter)>=StopNum)
                        break;
                }
                
                ShiftOrNot=1;
                NeedRoOrNOt=0;
            }
            
            if(RoTimes%RoTimes_PerShulffle==0){
                
                Parallel_Shuffle(CandiSeq,CandiNum,GRank0_StoreCandiSeq,TempVec,TempVec2,GRank,GroupSize,GWorld,Mode_ParaTech,InGRank,InGWorld);
                
                NeedRoOrNOt=0;
            }
            
            if(NeedRoOrNOt==1){
                
                Route2RouteVec(GRoute,RouteVec,CabsNumPerCore);
                Parallel_rotation(RouteVec,RouteVecTemp,GRank,GroupSize,Mode_ISend,GWorld);
                RouteVec2Route(GRoute,RouteVec,CabsNumPerCore);
                
                ShiftOrNot=1;
            }
            
            if(ShiftOrNot==1){
                
                ShiftTimes--;
                if(GRank+ShiftTimes<0)
                    ShiftTimes+=GroupSize;
                CurRouteIndex=CabsNumPerCore*((GRank+ShiftTimes)%GroupSize);
                GRoute_GDInf_GCabsPos_update(GRoute,GCabsPos,GDInf,Route,CabsPos,DInf,CurRouteIndex,CabsNumPerCore);
                
            }
            
            PTD_E0=PTD_MulX(GRoute,GDInf,SpotsPos,GCabsPos,CabsNumPerCore,DistVec,RouteProb);
            if(Mode_ParaTech==2)
                GRouteCandi_OriSoluVec_Ini(GRouteCandi,OriSoluVec,GRoute,CabsNumPerCore);
            
            // PTD_E0, GRoute, GRouteCandi, OriSoluVec, CandiSeq need to update again and again.
        }
    }
    
    *PTD_E0_Add=PTD_E0_AllRoute;
    *PTD_min_Add=PTD_min;
    *StepNum_Add=StepNum;
    
    for(i=0;i<CabsNumPerCore;i++){
        free(GRouteCandi_new[i]);
    }
    free(GRouteCandi_new);
    for(i=0;i<CabsNumPerCore;i++){
        free(GRouteCandi[i]);
    }
    free(GRouteCandi);
    free(ReplaceInfoVec);
    free(XVec2);
    free(XVec1);
    free(OriSoluVec);
    free(RouteVecTemp);
    free(RouteVec);
    for(i=0;i<CabsNumPerCore;i++){
        free(GCabsPos[i]);
    }
    free(GCabsPos);
    free(GDInf);
    for(i=0;i<CabsNumPerCore;i++)
        free(GRoute[i]);
    free(GRoute);
    for(i=0;i<CabsNumPerCore;i++)
        free(GRouteNew[i]);
    free(GRouteNew);
    free(TempVec2);
    free(TempVec);
    free(CandiSeq);
    for(i=0;i<CabsNum;i++){
        free(RouteNew[i]);
    }
    free(RouteNew);
    if(GRank==0)
        free(GRank0_StoreCandiSeq);
}



//Function Series Final: Overall Functions
void SAMethods(double** CabsPos, int** Route, const double* DInf, double** SpotsPos, const int SpotsNum, const int rank, const int size){
    int Mode_ParaTech; // 0 is no parallel, 1 is my techniques, 2 is my techniques + my MPeriod and MPattern, 3 is CDR, 4 is Lou.
    int Mode_ISend,Mode_PrintProc,Mode_CheckRoute; // 0 is no and 1 is yes.
    int ParaRecord,ParaRecordMax,RepeatedExpIndex,RepeatedExpInOneGroup;
    double TimeBeg,CompTime,timeTemp;
    int PrintPerSteps,RandomSeed;
    int StopNum; // Sequential parameters.
    double Temperature,Alpha; // Sequential parameters.
    int Steps_PerRo,RoTimes_PerInter,RoTimes_PerShulffle,SASteps_InInter; // My techniques parameters.
    int SuitRangeSteps; // My MPs
    double InRangeRate; // My MPs
    int MPeriod; // CDR or Lou's initial
    int MaxMP; // Lou
    int GRank,LGRank,InGRank,GroupSize,LGroupSize,InGroupSize; // My MPs. For 10 Cabs 2000 size, if LGroupSize is 50 and InGroupSize is 10, then 50 cores for 1 exp, 5=50/10 InGWorld each dealing with 2=10/5 cabs, and 10 cores dealing with 2 cabs. GWorld are across the InGWorld and GSize is 5. There are 10 GWorld.
    int StepNum,AccPeriod,OptRank,PrintFlag,i,j,k,*RouteVec,DupNum,DividNum;
    double PTD_E0,PTD_min,AccRateIni,PTD_ECheck,*DistVec,*RouteProb;
    FILE *ProcFile,*FinalFile;
    MPI_Comm GWorld,LGWorld,InGWorld;
    
    // ##### =>
    // SA Parameters:
    
    Mode_ParaTech=4;
    DividNum=1;
    
    Mode_CheckRoute=0;
    if(Mode_ParaTech==0){
        
        Mode_PrintProc=0;
        
        ParaRecordMax=1;
        RepeatedExpInOneGroup=1;
        GroupSize=1;
        
    }else if(Mode_ParaTech==1){
        
        Mode_ISend=1;
        Mode_PrintProc=0;
        
        ParaRecordMax=1;
        RepeatedExpInOneGroup=10;
        GroupSize=CabsNum/DividNum;
        
    }else if(Mode_ParaTech==2){
        
        Mode_ISend=1;
        Mode_PrintProc=0;
        
        ParaRecordMax=1;
        RepeatedExpInOneGroup=10;
        GroupSize=CabsNum; // GroupSize * InGroupSize is the LGroupSize. LGroupSize in PMode2 is the GroupSize in PMode1. GWorld here is each group that goes across the InGWorld.
        InGroupSize=2;
        
    }else if(Mode_ParaTech==3){
        
        Mode_PrintProc=0;
        
        ParaRecordMax=1;
        RepeatedExpInOneGroup=5;
        GroupSize=CabsNum/DividNum;
        
    }else if(Mode_ParaTech==4){
        
        Mode_PrintProc=0;
        
        ParaRecordMax=1;
        RepeatedExpInOneGroup=10;
        GroupSize=CabsNum/DividNum;
        
    }
    
    // <= #####
    
    if(Mode_ParaTech==2){
        
        // In Mode_ParaTech == 2, GRank, GSize, and GWorld are across the group. InGRank, InGSize, and InGWorld are inside the group. LGRank, LGSize, and LGWorld are the larger group.
        LGroupSize=GroupSize*InGroupSize;
        
        MPI_Comm_split(MPI_COMM_WORLD,rank/LGroupSize,0,&LGWorld);
        MPI_Comm_rank(LGWorld,&LGRank);
        
        MPI_Comm_split(LGWorld,LGRank/InGroupSize,0,&InGWorld);
        MPI_Comm_rank(InGWorld,&InGRank);
        
        MPI_Comm_split(LGWorld,LGRank%InGroupSize,0,&GWorld);
        MPI_Comm_rank(GWorld,&GRank);
        
    }else{
        MPI_Comm_split(MPI_COMM_WORLD,rank/GroupSize,0,&GWorld);
        MPI_Comm_rank(GWorld,&GRank);
    }
    
    timeTemp=MPI_Wtime();
    
    // ##### =>
    RandomSeed=(timeTemp-(int)timeTemp)*1000*(rank+1);//LIRed
    //RandomSeed=((10000*timeTemp-(int)(10000*timeTemp))*100)*(rank+1);//JN Intel
    //RandomSeed=(getpid()+1)*(MPI_Wtime()+0.0001)*10000*(rank+1); //SW
    // <= #####
    
    srand(RandomSeed);
    printf("Random seed: %d, rank: %d, timeTemp: %f.\n",RandomSeed,rank,timeTemp);
    fflush(stdout);
    
    if(Mode_PrintProc==1)
        ProcFile=fopen("proc.txt","a");
    FinalFile=fopen("FinalFile.txt","a");
    
    for(ParaRecord=0;ParaRecord<ParaRecordMax;ParaRecord++){
        for(RepeatedExpIndex=0;RepeatedExpIndex<RepeatedExpInOneGroup;RepeatedExpIndex++){
            TimeBeg=MPI_Wtime();
            
            // ##### =>
            // SA Parameters 2:
            
            PrintPerSteps=100000; // 100000 in MSR paper. Since Step%PrintPerSteps==1 printf, PrintPerSteps cannot be 1.
            StopNum=1000000; // 1000000 in MSR paper.
            
            if(Mode_ParaTech==0||Mode_ParaTech==3||Mode_ParaTech==4){
                
                Temperature=10.0; // 10 in MSR paper.
                Alpha=pow(1.0-pow(10,-6),sqrt(1.0/(double)(CabsNum)));
                //Alpha=pow(1.0-pow(10,-6),(1.0/(double)(ParaRecord+1)));
                // 1.0-pow(10,-6) for SA global and 1.0-pow(10,-7) for SA local in MSR paper.
                
                if(Mode_ParaTech==3){
                    
                    MPeriod=5; // Fixed MP in CDR.
                    
                }else if(Mode_ParaTech==4){
                    
                    MPeriod=5; // Initial MP in Lou.
                    MaxMP=1000; // Max MP in Lou.
                    
                }
                
            }else if(Mode_ParaTech==1||Mode_ParaTech==2){
                
                Temperature=10.0; // 10 in MSR paper.
                Alpha=pow(1.0-pow(10,-6),sqrt(1.0/(double)(CabsNum)));
                // 1.0-pow(10,-6) for SA global and 1.0-pow(10,-7) for SA local in MSR paper.
                
                if(Mode_ParaTech==1){
                    
                    Steps_PerRo=pow(10,1+0); // 1 + 0
                    SASteps_InInter=pow(10,1+0); // should be 0 + 0. But want it to be 0 + 1.
                    RoTimes_PerInter=CabsNum*(0+1); // 0 + 1
                    RoTimes_PerShulffle=CabsNum*(0+1); // 0 + 1
                    
                }else if(Mode_ParaTech==2){
                    
                    Steps_PerRo=10; // 10
                    SASteps_InInter=10;
                    RoTimes_PerInter=CabsNum;
                    RoTimes_PerShulffle=CabsNum;
                    
                    SuitRangeSteps=pow(10,1); // Named SearchSteps in the functions.
                    
                    InRangeRate=1.0/(double)(InGroupSize-1); // May need to change for the final version.
                    MPeriod=5;
                    AccPeriod=1000;
                    AccRateIni=1.0;
                }
                
            }
            
            // <= #####
            
            if(Mode_ParaTech==0){
                
                Seq_only(CabsPos,Route,DInf,SpotsPos,SpotsNum,Mode_PrintProc,PrintPerSteps,StopNum,Temperature,Alpha,ProcFile,&PTD_E0,&PTD_min,&StepNum);
                // Output: Route, ProcFile, PTD_E0, PTD_min, StepNum.
                
            }else if(Mode_ParaTech==3||Mode_ParaTech==4){
                
                CDR_Lou(CabsPos,Route,DInf,SpotsPos,SpotsNum,Mode_PrintProc,PrintPerSteps,StopNum,Temperature,Alpha,ProcFile,&PTD_E0,&PTD_min,&StepNum,Mode_ParaTech,MPeriod,MaxMP,GRank,GroupSize,GWorld);
                // Output: Route, ProcFile, PTD_E0, PTD_min, StepNum.
                
            }else if(Mode_ParaTech==1||Mode_ParaTech==2){
                
                MyTech_MyMPs(CabsPos,Route,DInf,SpotsPos,SpotsNum,Mode_PrintProc,PrintPerSteps,StopNum,Temperature,Alpha,ProcFile,&PTD_E0,&PTD_min,&StepNum,Mode_ParaTech,GRank,GroupSize,GWorld,LGRank,LGroupSize,LGWorld,InGRank,InGroupSize,InGWorld,Mode_ISend,Steps_PerRo,RoTimes_PerInter,RoTimes_PerShulffle,SASteps_InInter,SuitRangeSteps,InRangeRate,MPeriod,AccPeriod,AccRateIni,&OptRank);
                // Output: OptRank (Mode_ParaTech1: GRank, Mode_ParaTech2: LGRank), OptRank's Route (that has PTD_E0), ProcFile, PTD_E0, PTD_min, StepNum.
            }
            
            CompTime=MPI_Wtime()-TimeBeg;
            
            PrintFlag=0;
            if(Mode_ParaTech==0){
                PrintFlag=1;
            }else if(Mode_ParaTech==3||Mode_ParaTech==4){
                if(GRank==0)
                    PrintFlag=1;
            }else if(Mode_ParaTech==1){
                if(GRank==OptRank)
                    PrintFlag=1;
            }else if(Mode_ParaTech==2){
                if(LGRank==OptRank)
                    PrintFlag=1;
            }
            
            if(PrintFlag==1){
                
                fprintf(FinalFile,"%f %f %d %f %d %d %d %d\n",PTD_E0,PTD_min,StepNum,CompTime,Mode_ParaTech,ParaRecord,RepeatedExpIndex,RandomSeed);
                fflush(FinalFile);
                
                if(Mode_CheckRoute==1){
                    
                    RouteVec=(int*)malloc(CabsNum*RouteLimit*sizeof(int));
                    DistVec=(double*)malloc(RouteLimit*sizeof(double));
                    RouteProb=(double*)malloc(RouteLimit*sizeof(double));
                    
                    PTD_ECheck=PTD_MulX(Route,DInf,SpotsPos,CabsPos,CabsNum,DistVec,RouteProb);
                    
                    printf("Route:\n");
                    for(i=0,k=0;i<CabsNum;i++){
                        for(j=1;j<=Route[0][0];j++){
                            printf("%d ",Route[i][j]);
                        }
                        printf("\n");
                    }
                    printf("%f %f\n",PTD_E0,PTD_ECheck);
                    
                    Route2RouteVec(Route,RouteVec,CabsNum);
                    SortInt_NoPlace(RouteVec,CabsNum*RouteLimit);
                    
                    for(i=1,DupNum=0;i<CabsNum*RouteLimit;i++){
                        if(RouteVec[i]==RouteVec[i-1])
                            DupNum++;
                    }
                    if(DupNum!=CabsNum-1){
                        if(DupNum==CabsNum)
                            printf("Maybe Error! Maybe not because Pos NO. %d may also be in the route\n",Route[0][0]);
                        else
                            printf("Error! The route results are unexpected.\n");
                    }
                        
                    
                    for(i=0;i<CabsNum*RouteLimit;i++)
                        printf("%d ",RouteVec[i]);
                    printf("\n");
                    
                    free(RouteProb);
                    free(DistVec);
                    free(RouteVec);
                }
                
                // PTD_E0's Route is here.
                // RandomSeed is appropriate for repeating experiments but we need to know which rank it is from and WallTime should not be that sensitive. For example, precision to second can be more stable. If so, remember to truncate the time below second. Or 1.273873 * 1000 is unpredicable.
            }
        }
    }
    
    fclose(FinalFile);
    if(Mode_PrintProc==1)
        fclose(ProcFile);
    
}

double IBPMethod(const double* StartPos, int* OptRoute, const double Dmax, double** PPPosProb, const int PPNum, const int RouteToChoose){
    
    //Algorithm 1:
    int i,j,k,RouteL,RouteNo,PPNo,flag_duplicate,RL_routeNum_temp,flag_check,NotEmpty,RouteN,RouteNo_check,SourcePtTemp;
    int *r_vec,*q_vec,***RL_route_mat,*RL_routeNum;
    double F1_r_vec,PE_r_vec,F1_rpre_vec,PE_rpre_vec;
    double ***RL_FP_mat,**DistMat;
    double TimeBegin,TimeBegin2,TimeEnd;
    //Algorithm 2:
    int **RCLIndex_mat,*RCLIndex_vec,*RouteEnd,*RouteIni,souce_c,IsGoodRoute,RouteNo_RCL,r_index,q_index,**BPRL_route_mat,SelectedRouteNum,Index_RL,count;
    double F1_r,PE_r,F1_q,PE_q,**BPRL_FP_mat;
    
    //Algorithm 3:
    int OpRL_route_index;
    double Fmin,Ftemp,DistC0C;
    
    double OptVal[20];
    
        
    RouteN=RouteToChoose;
    
    r_vec=(int*)malloc(RouteN*sizeof(int));
    q_vec=(int*)malloc(RouteN*sizeof(int));
    RL_routeNum=(int*)malloc(RouteN*sizeof(int));
    
    //Assign initial value of vector RL_routeNum
    RL_routeNum[0]=PPNum;
    for(i=1;i<RouteN;i++)
    {
        RL_routeNum[i]=RL_routeNum[i-1]*(PPNum-i)/i;
    }
    
    RL_route_mat=(int***)malloc(RouteN*sizeof(int**));
    for(i=0;i<RouteN;i++)
    {
        RL_route_mat[i]=(int**)malloc(RL_routeNum[i]*sizeof(int*));
        for(j=0;j<RL_routeNum[i];j++)
        {
            RL_route_mat[i][j]=(int*)malloc((i+1)*sizeof(int));
        }
    }
    
    RL_FP_mat=(double***)malloc(RouteN*sizeof(double**));
    for(i=0;i<RouteN;i++)
    {
        RL_FP_mat[i]=(double**)malloc(RL_routeNum[i]*sizeof(double*));
        for(j=0;j<RL_routeNum[i];j++)
        {
            RL_FP_mat[i][j]=(double*)malloc(2*sizeof(double));
        }
    }
    //Assign initial value of RL_FP_mat and RL_route_mat
    for(i=0;i<PPNum;i++)
    {
        RL_route_mat[0][i][0]=i;
        RL_FP_mat[0][i][0]=0.0;
        RL_FP_mat[0][i][1]=PPPosProb[i][2];
    }
    
    //Allocate and assign values of DistMat
    DistMat=(double**)malloc(PPNum*sizeof(double*));
    for(i=0;i<PPNum;i++)
    {
        DistMat[i]=(double*)malloc(PPNum*sizeof(double));
    }
    for(i=0;i<PPNum;i++)
    {
        for(j=0;j<PPNum;j++)
        {
            DistMat[i][j]=sqrt((PPPosProb[i][0]-PPPosProb[j][0])*(PPPosProb[i][0]-PPPosProb[j][0])+(PPPosProb[i][1]-PPPosProb[j][1])*(PPPosProb[i][1]-PPPosProb[j][1]));
        }
    }
    
    //Main iteration starts
    for(RouteL=1;RouteL<RouteN;RouteL++)
    {
        RL_routeNum_temp=0;
        //For each route in previous optimal set
        for(RouteNo=0;RouteNo<RL_routeNum[RouteL-1];RouteNo++)
        {
            //For each pick-up point
            for(PPNo=0;PPNo<PPNum;PPNo++)
            {
                //Check if it belongs to the pick up points in that sepcfic route
                flag_duplicate=0;
                for(j=0;j<RouteL;j++)
                {
                    if(PPNo==RL_route_mat[RouteL-1][RouteNo][j])
                    {
                        flag_duplicate=1;
                        break;
                    }
                }
                if(flag_duplicate==1)
                {
                    continue;
                }
                //If it is here,flag_duplicate==0, indicating PPNo is not the same as any of the points there
                //Then we can form a new route stored in the last several positions of r_vec
                r_vec[RouteN-RouteL-1]=PPNo;
                for(j=RouteN-RouteL;j<RouteN;j++)
                {
                    r_vec[j]=RL_route_mat[RouteL-1][RouteNo][j-(RouteN-RouteL)];
                }
                //Find source point
                SourcePtTemp=r_vec[RouteN-RouteL];
                //Use the info before and some info from the current pppoint to compute the F1 and PE
                F1_rpre_vec=RL_FP_mat[RouteL-1][RouteNo][0];
                PE_rpre_vec=RL_FP_mat[RouteL-1][RouteNo][1];
                F1_r_vec=F1_rpre_vec*(1-PPPosProb[SourcePtTemp][2])+DistMat[PPNo][SourcePtTemp]*PE_rpre_vec;
                PE_r_vec=PE_rpre_vec*(1-PPPosProb[PPNo][2])+PPPosProb[PPNo][2];
                //Show the optimal route recursively
                NotEmpty=0;
                for(RouteNo_check=0;RouteNo_check<RL_routeNum_temp;RouteNo_check++)
                {
                    //Check if the sourvce point of the new route is the same as the one in an optimal route
                    if(r_vec[RouteN-RouteL-1]==RL_route_mat[RouteL][RouteNo_check][0])
                    {
                        //Check if the two routes have the same pick-up points
                        for(j=RouteN-RouteL;j<RouteN;j++)
                        {
                            flag_check=0;
                            for(k=1;k<RouteL+1;k++)
                            {
                                if(r_vec[j]==RL_route_mat[RouteL][RouteNo_check][k])
                                {
                                    flag_check++;
                                    break;
                                }
                            }
                            if(flag_check==0)
                                break;
                        }
                        if(flag_check==0)
                            continue;
                        //If it is here, the two routes have the same pick-ip points
                        NotEmpty++;
                        //If the new route is better than the optimal route, replace it
                        if(F1_r_vec<RL_FP_mat[RouteL][RouteNo_check][0])
                        {
                            RL_FP_mat[RouteL][RouteNo_check][0]=F1_r_vec;
                            RL_FP_mat[RouteL][RouteNo_check][1]=PE_r_vec;
                            for(j=0;j<RouteL+1;j++)
                            {
                                RL_route_mat[RouteL][RouteNo_check][j]=r_vec[j+RouteN-RouteL-1];
                            }
                            //break;
                            //If the else below does not happen, a "break" should be added.
                        }else if(F1_r_vec==RL_FP_mat[RouteL][RouteNo_check][0])
                        {
                            //printf("!!Error. Two values are the same.");//?????
                        }
                        //break;
                        //At most one route in the optimal route set has the same source and pick-up points as the new route. After we examine this route, there will be no other routes satisfying this condition. And we don't need to check the other optimal routes.
                    }
                }
                //If NotEmpty==0, there is no optimal route having the same source and pick-up points as the new route. Then this route CANNOT BE ILIMINATED. We keep it.
                if(NotEmpty==0)
                {
                    RL_FP_mat[RouteL][RL_routeNum_temp][0]=F1_r_vec;
                    RL_FP_mat[RouteL][RL_routeNum_temp][1]=PE_r_vec;
                    for(j=0;j<RouteL+1;j++)
                    {
                        RL_route_mat[RouteL][RL_routeNum_temp][j]=r_vec[j+RouteN-RouteL-1];
                    }
                    RL_routeNum_temp++;
                }
                
            }
        }
    }

    //Algorithm 1 completes.
    //Algorithm 2 starts.
    
    for(RouteL=0;RouteL<RouteToChoose;RouteL++){
        //Assign values for RouteL !!!
        //RouteL=RouteToChoose-1;//a-1, a is the length.
        
        
        RCLIndex_mat=(int**)malloc(PPNum*sizeof(int*));
        for(i=0;i<PPNum;i++)
        {
            RCLIndex_mat[i]=(int*)malloc(RL_routeNum[RouteL]/PPNum*sizeof(int));
        }
        RCLIndex_vec=(int*)malloc(PPNum*sizeof(int));
        RouteIni=(int*)malloc(PPNum*sizeof(int));
        RouteEnd=(int*)malloc(PPNum*sizeof(int));
        
        for(i=0;i<PPNum;i++)
        {
            RCLIndex_vec[i]=0;
        }
        for(i=0;i<RL_routeNum[RouteL];i++)
        {
            souce_c=RL_route_mat[RouteL][i][0];
            RCLIndex_mat[souce_c][RCLIndex_vec[souce_c]++]=i;
        }
        for(PPNo=0;PPNo<PPNum;PPNo++)
        {
            RouteEnd[PPNo]=RL_routeNum[RouteL]/PPNum;
            RouteIni[PPNo]=0;
            for(RouteNo=1;RouteNo<RouteEnd[PPNo];RouteNo++)
            {
                IsGoodRoute=0;
                for(RouteNo_RCL=RouteNo-1;RouteNo_RCL>=RouteIni[PPNo];RouteNo_RCL--)
                {
                    //New solution is good
                    r_index=RCLIndex_mat[PPNo][RouteNo];
                    q_index=RCLIndex_mat[PPNo][RouteNo_RCL];
                    F1_r=RL_FP_mat[RouteL][r_index][0];
                    PE_r=RL_FP_mat[RouteL][r_index][1];
                    F1_q=RL_FP_mat[RouteL][q_index][0];
                    PE_q=RL_FP_mat[RouteL][q_index][1];
                    if((F1_r<F1_q&&PE_q<=PE_r)||(F1_r<=F1_q&&PE_q<PE_r))
                    {
                        IsGoodRoute=1;
                        RCLIndex_mat[PPNo][RouteNo_RCL]=RCLIndex_mat[PPNo][RouteIni[PPNo]];
                        RouteIni[PPNo]++;
                        RouteNo_RCL++;
                    }else if(IsGoodRoute==0&&((F1_q<F1_r&&PE_r<=PE_q)||(F1_q<=F1_r&&PE_r<PE_q)))//New solution is bad
                    {
                        RCLIndex_mat[PPNo][RouteNo]=RCLIndex_mat[PPNo][RouteEnd[PPNo]-1];
                        RouteEnd[PPNo]--;
                        RouteNo--;
                        break;
                    }
                }
            }
        }
        
        SelectedRouteNum=0;
        for(i=0;i<PPNum;i++){
            SelectedRouteNum+=RouteEnd[i]-RouteIni[i];
        }
        BPRL_route_mat=(int**)malloc(SelectedRouteNum*sizeof(int*));
        for(i=0;i<SelectedRouteNum;i++){
            BPRL_route_mat[i]=(int*)malloc((RouteL+1)*sizeof(int));
        }
        BPRL_FP_mat=(double**)malloc(SelectedRouteNum*sizeof(double*));
        for(i=0;i<SelectedRouteNum;i++){
            BPRL_FP_mat[i]=(double*)malloc(2*sizeof(double));
        }
        
        count=0;
        for(i=0;i<PPNum;i++){
            for(j=RouteIni[i];j<RouteEnd[i];j++){
                Index_RL=RCLIndex_mat[i][j];
                BPRL_FP_mat[count][0]=RL_FP_mat[RouteL][Index_RL][0];
                BPRL_FP_mat[count][1]=RL_FP_mat[RouteL][Index_RL][1];
                for(k=0;k<RouteL+1;k++){
                    BPRL_route_mat[count][k]=RL_route_mat[RouteL][Index_RL][k];
                }
                count++;
            }
        }
        //Algorithm 2 ends
        //Algorithm 3 starts
        
        OpRL_route_index=0;
        SourcePtTemp=BPRL_route_mat[0][0];
        DistC0C=sqrt((PPPosProb[SourcePtTemp][0]-StartPos[0])*(PPPosProb[SourcePtTemp][0]-StartPos[0])+(PPPosProb[SourcePtTemp][1]-StartPos[1])*(PPPosProb[SourcePtTemp][1]-StartPos[1]));
        Fmin=BPRL_FP_mat[0][0]*(1.0-PPPosProb[SourcePtTemp][2])+DistC0C*BPRL_FP_mat[0][1]+Dmax*(1-BPRL_FP_mat[0][1]);
        
        
        for(RouteNo=1;RouteNo<SelectedRouteNum;RouteNo++)
        {
            SourcePtTemp=BPRL_route_mat[RouteNo][0];
            DistC0C=sqrt((PPPosProb[SourcePtTemp][0]-StartPos[0])*(PPPosProb[SourcePtTemp][0]-StartPos[0])+(PPPosProb[SourcePtTemp][1]-StartPos[1])*(PPPosProb[SourcePtTemp][1]-StartPos[1]));
            Ftemp=BPRL_FP_mat[RouteNo][0]*(1.0-PPPosProb[SourcePtTemp][2])+DistC0C*BPRL_FP_mat[RouteNo][1]+Dmax*(1-BPRL_FP_mat[RouteNo][1]);
            if(Ftemp<Fmin){
                Fmin=Ftemp;
                OpRL_route_index=RouteNo;
            }
        }

        OptVal[RouteL]=Fmin;
        //ALgorithm 3 ends
        //Algorithm 2 starts
        
        
        for(i=0;i<PPNum;i++)
        {
            free(RCLIndex_mat[i]);
        }
        free(RCLIndex_mat);
        free(RCLIndex_vec);
        free(RouteIni);
        free(RouteEnd);
        
        for(i=0;i<SelectedRouteNum;i++){
            free(BPRL_route_mat[i]);
        }
        free(BPRL_route_mat);
        
        for(i=0;i<SelectedRouteNum;i++){
            free(BPRL_FP_mat[i]);
        }
        free(BPRL_FP_mat);
        //Algorithm 2 ends.
    }
    
    //Algorithm 1 starts.
    
    
    free(r_vec);
    free(q_vec);
    
    for(i=0;i<RouteN;i++)
    {
        for(j=0;j<RL_routeNum[i];j++)
        {
            free(RL_route_mat[i][j]);
        }
        free(RL_route_mat[i]);
    }
    free(RL_route_mat);
    
    for(i=0;i<RouteN;i++)
    {
        for(j=0;j<RL_routeNum[i];j++)
        {
            free(RL_FP_mat[i][j]);
        }
        free(RL_FP_mat[i]);
    }
    free(RL_FP_mat);
    
    free(RL_routeNum);
    for(i=0;i<PPNum;i++)
    {
        free(DistMat[i]);
    }
    free(DistMat);
    
    return OptVal[RouteToChoose-1];
    
}

double EnuMethod(const double* IniPos, int* OptRoute, const double Dmax, double** ppPos, const int ppNum, const int RouteL){
    
    int i,j,iter_m,flag_temp,flag_toBreak,PTD_CompNum;
    int *CurRoute;
    double CurrentPTD,MinPTD;
    double *CurPTD,*CurProb,*DistVec;
    double **DistMAT;
    
    DistMAT=(double**)malloc(ppNum*sizeof(double*));
    for(i=0;i<ppNum;i++)
    {
        DistMAT[i]=(double*)malloc(ppNum*sizeof(double));
    }
    for(i=0;i<ppNum;i++)
    {
        for(j=0;j<ppNum;j++)
        {
            DistMAT[i][j]=sqrt((ppPos[i][0]-ppPos[j][0])*(ppPos[i][0]-ppPos[j][0])+(ppPos[i][1]-ppPos[j][1])*(ppPos[i][1]-ppPos[j][1]));
        }
    }
    
    CurPTD=(double*)malloc(RouteL*sizeof(double));
    CurProb=(double*)malloc(RouteL*sizeof(double));
    DistVec=(double*)malloc(RouteL*sizeof(double));
    CurRoute=(int*)malloc(RouteL*sizeof(int));
    
    MinPTD=10000000000.0;
    
    for(i=0;i<RouteL;i++){
        CurRoute[i]=-1;
    }
    iter_m=0;
    flag_toBreak=0;
    PTD_CompNum=0;
    
    while(1){
        
        while(1){
            
            while(1){
                CurRoute[iter_m]++;
                if(CurRoute[iter_m]==ppNum){
                    if(iter_m==0){
                        flag_toBreak=1;
                        break;
                    }
                    CurRoute[iter_m]=-1;
                    iter_m--;
                }else{
                    break;
                }
            }
            if(flag_toBreak==1)
                break;
            
            
            flag_temp=0;
            for(i=0;i<iter_m;i++){
                if(CurRoute[iter_m]==CurRoute[i]){
                    flag_temp=1;
                    break;
                }
            }
            if(flag_temp==0)
                break;
        }
        if(flag_toBreak==1)
            break;
        
        if(iter_m==0){
            CurProb[iter_m]=1-ppPos[CurRoute[iter_m]][2];
            DistVec[iter_m]=sqrt((IniPos[0]-ppPos[CurRoute[i]][0])*(IniPos[0]-ppPos[CurRoute[i]][0])+(IniPos[1]-ppPos[CurRoute[i]][1])*(IniPos[1]-ppPos[CurRoute[i]][1]));
            CurPTD[iter_m]=DistVec[iter_m]*ppPos[CurRoute[iter_m]][2];
        }else{
            CurProb[iter_m]=CurProb[iter_m-1]*(1-ppPos[CurRoute[iter_m]][2]);
            DistVec[iter_m]=DistVec[iter_m-1]+DistMAT[CurRoute[iter_m-1]][CurRoute[iter_m]];
            CurPTD[iter_m]=CurPTD[iter_m-1]+CurProb[iter_m-1]*ppPos[CurRoute[iter_m]][2]*DistVec[iter_m];
        }
        
        if(iter_m==RouteL-1){
            PTD_CompNum++;
            CurrentPTD=CurPTD[iter_m]+CurProb[iter_m]*Dmax;
            if(CurrentPTD<MinPTD){
                MinPTD=CurrentPTD;
                for(i=0;i<RouteL;i++){
                    OptRoute[i]=CurRoute[i];
                }
            }
            
        }else{
            iter_m++;
        }
    }
    
    for(i=0;i<ppNum;i++)
        free(DistMAT[i]);
    free(DistMAT);
    free(CurPTD);
    free(CurProb);
    free(DistVec);
    free(CurRoute);
    
    return MinPTD;
    
}

void IBP_Enu(double** CabsPos, int** Route, const double* DInf, double** SpotsPos, const int SpotsNum, const int rank, const int size, const int Mode_Method){
    int i,j,*Place,**PlaceAllCores,*PlaceAllCores_Vec,**Subset,*SubsetVec,SubsetSize,flag,temp,*Record,RecordSize,i0,*LocalSubset,*SingleRoute,*MyRoute,*Index;
    double *DistToStart,**ppPos,LocalPTD,AllPTD,TimeBeg,CompTime,tempVal;
    
    // ##### =>
    SubsetSize=17;
    // <= #####
    
    DistToStart=(double*)malloc(SpotsNum*sizeof(double));
    Place=(int*)malloc(SpotsNum*sizeof(int));
    LocalSubset=(int*)malloc(SubsetSize*sizeof(int));
    SingleRoute=(int*)malloc(Route[0][0]*sizeof(int));
    MyRoute=(int*)malloc(Route[0][0]*sizeof(int));
    
    ppPos=(double**)malloc(SubsetSize*sizeof(double*));
    for(i=0;i<SubsetSize;i++){
        ppPos[i]=(double*)malloc(3*sizeof(double));
    }
    
    if(rank==0){
        PlaceAllCores=(int**)malloc(CabsNum*sizeof(int*));
        for(i=0;i<CabsNum;i++){
            PlaceAllCores[i]=(int*)malloc(CabsNum*SubsetSize*sizeof(int));
        }
        PlaceAllCores_Vec=(int*)malloc(CabsNum*CabsNum*SubsetSize*sizeof(int));
        Subset=(int**)malloc(CabsNum*sizeof(int*));
        for(i=0;i<CabsNum;i++){
            Subset[i]=(int*)malloc(SubsetSize*sizeof(int));
        }
        SubsetVec=(int*)malloc(CabsNum*SubsetSize*sizeof(int));
        Index=(int*)malloc(CabsNum*sizeof(int));
        Record=(int*)malloc(CabsNum*SubsetSize*sizeof(int));
    }
    
    TimeBeg=MPI_Wtime();
    
    for(i=0;i<SpotsNum;i++){
        DistToStart[i]=pow(pow(SpotsPos[i][0]-CabsPos[rank][0],2)+pow(SpotsPos[i][1]-CabsPos[rank][1],2),0.5);
        Place[i]=i;
    }
    
    for(i=0;i<CabsNum*SubsetSize;i++){
        for(j=SpotsNum-1;j>i;j--){
            if(DistToStart[j-1]>DistToStart[j]){
                tempVal=DistToStart[j-1];
                DistToStart[j-1]=DistToStart[j];
                DistToStart[j]=tempVal;
                
                SwapInt(&Place[j-1],&Place[j]);
            }
        }
    }
    
    MPI_Gather(Place,CabsNum*SubsetSize,MPI_INT,PlaceAllCores_Vec,CabsNum*SubsetSize,MPI_INT,0,MPI_COMM_WORLD);
    
    if(rank==0){
        for(i=0;i<CabsNum;i++){
            for(j=0;j<CabsNum*SubsetSize;j++){
                PlaceAllCores[i][j]=PlaceAllCores_Vec[i*CabsNum*SubsetSize+j];
            }
        }
        
        for(i=0;i<CabsNum;i++)
            Index[i]=0;
        
        RecordSize=0;
        for(j=0;j<SubsetSize;j++){
            for(i=0;i<CabsNum;i++){
                flag=1;
                while(flag){
                    temp=PlaceAllCores[i][Index[i]++];
                    flag=0;
                    for(i0=0;i0<RecordSize;i0++){
                        if(Record[i0]==temp){
                            flag=1;
                            break;
                        }
                    }
                }
                
                Subset[i][j]=temp;
                Record[RecordSize++]=temp;
            }
        }
        
        for(i=0;i<CabsNum;i++){
            for(j=0;j<SubsetSize;j++){
                SubsetVec[i*SubsetSize+j]=Subset[i][j];
            }
        }
        
    }
    
    MPI_Scatter(SubsetVec,SubsetSize,MPI_INT,LocalSubset,SubsetSize,MPI_INT,0,MPI_COMM_WORLD);
    
    for(i=0;i<SubsetSize;i++){
        ppPos[i][0]=SpotsPos[LocalSubset[i]][0];
        ppPos[i][1]=SpotsPos[LocalSubset[i]][1];
        ppPos[i][2]=SpotsPos[LocalSubset[i]][2];
    }
    
    if(Mode_Method==1)
        LocalPTD=IBPMethod(CabsPos[rank],SingleRoute,DInf[rank],ppPos,SubsetSize,Route[0][0]);
    else if(Mode_Method==2)
        LocalPTD=EnuMethod(CabsPos[rank],SingleRoute,DInf[rank],ppPos,SubsetSize,Route[0][0]);
    
    if(Mode_Method==2){
        for(i=0;i<Route[0][0];i++){
            MyRoute[i]=LocalSubset[SingleRoute[i]];
        }
    }
    
    MPI_Reduce(&LocalPTD,&AllPTD,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    
    CompTime=MPI_Wtime()-TimeBeg;
    
    if(rank==0){
        printf("%f %f\n",AllPTD,CompTime);
    }
    
    if(rank==0){
        free(Record);
        free(Index);
        free(SubsetVec);
        for(i=0;i<CabsNum;i++){
            free(Subset[i]);
        }
        free(Subset);
        free(PlaceAllCores_Vec);
        for(i=0;i<CabsNum;i++){
            free(PlaceAllCores[i]);
        }
        free(PlaceAllCores);
    }
    
    for(i=0;i<SubsetSize;i++)
        free(ppPos[i]);
    free(ppPos);
    free(MyRoute);
    free(SingleRoute);
    free(LocalSubset);
    free(Place);
    free(DistToStart);
}

void KDDPaper(const int rank, const int size){
    int i,j,k;
    double *PosRecord_temp,temp_input;
    
    int RouteLength,SpotsNum,**Route;
    double **CabsPos,DInf[CabsNum],**SpotsPos;
    
    int row,col,Mode_Method,Mode_Data;
    double DeltaR,DeltaC;
    
    FILE *SpotsPosFile;
    
    CabsPos=(double**)malloc(CabsNum*sizeof(double*));
    for(i=0;i<CabsNum;i++){
        CabsPos[i]=(double*)malloc(2*sizeof(double));
    }
    
    // ##### =>
    // OF Parameters:
    Mode_Method=0; // 0 is PSA, 1 is IBP, and 2 is Enu.      For IBP and Enu, size == CabsNum.
    Mode_Data=0; // 0 is RD_BJ, 1 is Syn.
    
    if(CabsNum==6){
        CabsPos[0][0]=0.50;
        CabsPos[0][1]=0.50;
        CabsPos[1][0]=0.51;
        CabsPos[1][1]=0.51;
        CabsPos[2][0]=0.52;
        CabsPos[2][1]=0.52;
        CabsPos[3][0]=0.53;
        CabsPos[3][1]=0.53;
        CabsPos[4][0]=0.54;
        CabsPos[4][1]=0.54;
        CabsPos[5][0]=0.55;
        CabsPos[5][1]=0.55;
        DInf[0]=1.0;
        DInf[1]=1.0;
        DInf[2]=1.0;
        DInf[3]=1.0;
        DInf[4]=1.0;
        DInf[5]=1.0;
    }else if(CabsNum==96){
        row=12;
        col=8;
        DeltaR=0.05/((double)row-1.0);
        DeltaC=0.05/((double)col-1.0);
        for(i=0,k=0;i<row;i++){
            for(j=0;j<col;j++){
                DInf[k]=1.0;
                CabsPos[k][0]=0.5+DeltaR*i;
                CabsPos[k++][1]=0.5+DeltaC*j;
            }
        }
        
        if(rank==0){
            for(i=0;i<CabsNum;i++){
                printf("%f ",CabsPos[i][0]);
                printf("%f ",CabsPos[i][1]);
                printf("%f\n",DInf[i]);
            }
        }
        
    }
    
    // <= #####
    
    RouteLength=RouteLimit-1;
    if(rank==0){
        if(Mode_Data==0){
            SpotsPosFile=fopen("RD_BJ_21824.txt","r");
            fscanf(SpotsPosFile,"%d",&SpotsNum);
        }else if(Mode_Data==1){
            SpotsPosFile=fopen("SynData.txt","r");
            fscanf(SpotsPosFile,"%d",&SpotsNum);
        }
    }
    
    MPI_Bcast(&SpotsNum,1,MPI_INT,0,MPI_COMM_WORLD);
    PosRecord_temp=(double*)malloc(3*SpotsNum*sizeof(double));
    SpotsPos=(double**)malloc(SpotsNum*sizeof(double*));
    for(i=0;i<SpotsNum;i++){
        SpotsPos[i]=(double*)malloc(3*sizeof(double));
    }
    Route=(int**)malloc(CabsNum*sizeof(int*));
    for(i=0;i<CabsNum;i++){
        Route[i]=(int*)malloc(RouteLimit*sizeof(int));
    }
    
    if(rank==0){
        for(i=0;i<3*SpotsNum;i++){
            fscanf(SpotsPosFile,"%lf",&temp_input);
            PosRecord_temp[i]=temp_input;
        }
        fclose(SpotsPosFile);
    }
    MPI_Bcast(PosRecord_temp,3*SpotsNum,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for(i=0;i<SpotsNum;i++){
        SpotsPos[i][0]=PosRecord_temp[3*i+0];
        SpotsPos[i][1]=PosRecord_temp[3*i+1];
        SpotsPos[i][2]=PosRecord_temp[3*i+2];
    }
    free(PosRecord_temp);
    for(i=0;i<CabsNum;i++){
        Route[i][0]=RouteLength;
    }
    
    if(Mode_Method==0)
        SAMethods(CabsPos,Route,DInf,SpotsPos,SpotsNum,rank,size);
    else
        IBP_Enu(CabsPos,Route,DInf,SpotsPos,SpotsNum,rank,size,Mode_Method);
    
    
    for(i=0;i<CabsNum;i++){
        free(Route[i]);
    }
    free(Route);
    for(i=0;i<SpotsNum;i++)
        free(SpotsPos[i]);
    free(SpotsPos);
    for(i=0;i<CabsNum;i++){
        free(CabsPos[i]);
    }
    free(CabsPos);
}

int main(int argc, char** argv)
{
    int rank,size;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    KDDPaper(rank,size);
    
    MPI_Finalize();
    return 0;
}


