#include<cstdio>
#include<vector>
#include<list>
#include<cmath>
#include<omp.h>
using namespace std;

#define THREADS 32
#define MAX 10000
#define NOT_CONNECTED -1

int nodesCount;
int FWDist[MAX][MAX];

struct Edge
{   
    int to,len;
};


int main()
{
    //double st,et,etr;
    //st = omp_get_wtime();
    #pragma omp parallel for num_threads(32)
    for (int i=0;i<MAX;++i){
        for (int j=0;j<MAX;++j){
            FWDist[i][j]=NOT_CONNECTED;

        }
        FWDist[i][i]=0;
    }

    scanf("%d", &nodesCount);
    vector<list<Edge> > adjlist(nodesCount,list<Edge>());

    int a, b, c, p, edgesCount = 0;
    while(scanf("%d %d %d", &a, &b, &c)!= EOF){
        if ( a > nodesCount || b > nodesCount){
            printf("Vertex index out of boundary.");
            return -1;
        }
        edgesCount++;
        FWDist[a][b]=c;

        Edge tmp;
        p = a-1;
        tmp.to = b-1;
        tmp.len = c;
        adjlist[p].push_back(tmp);

    }
    //etr = omp_get_wtime();
    //printf("read data time: %lf s\n",etr-st);
    int diameter = -1;
    if(edgesCount < 0.03*nodesCount*nodesCount){
        #pragma omp parallel for num_threads(32)
        for(int i=0; i<nodesCount; ++i){
            vector<int> dist,path;
            int beg = i;

            const int INF=0x7FFFFFFF,NODE=adjlist.size();
            dist.assign(NODE,INF);
            path.assign(NODE,-1);
            list<int> que(1,beg);
            vector<int> cnt(NODE,0);
            vector<bool> flag(NODE,0);
            dist[beg]=0;
            cnt[beg]=flag[beg]=1;
            while(!que.empty())
            {
                const int now=que.front();
                que.pop_front();
                flag[now]=0;
                for(list<Edge>::const_iterator
                        i=adjlist[now].begin(); i!=adjlist[now].end(); ++i)
                    if(dist[i->to]>dist[now]+i->len)
                    {
                        dist[i->to]=dist[now]+i->len;
                        path[i->to]=now;
                        if(!flag[i->to])
                        {
                            //if(NODE==++cnt[i->to])return 1;
                            if(!que.empty()&&dist[i->to]<dist[que.front()])
                                que.push_front(i->to);
                            else que.push_back(i->to);
                            flag[i->to]=1;
                        }
                    }
            }

            for(int j=0; j<nodesCount; ++j)
            {
                if(dist[j]<2000000000 && dist[j] > diameter)
                {
                    diameter = dist[j];
                }
            }
        }
    }
    else{
        #pragma omp parallel for collapse(2) num_threads(32)
        for (int k=1;k<=nodesCount;++k){
            for (int i=1;i<=nodesCount;++i){
                if (FWDist[i][k]!=NOT_CONNECTED){
                    for (int j=1;j<=nodesCount;++j){
                        if (FWDist[k][j]!=NOT_CONNECTED && (FWDist[i][j]==NOT_CONNECTED || FWDist[i][k]+FWDist[k][j]<FWDist[i][j])){
                            FWDist[i][j]=FWDist[i][k]+FWDist[k][j];
                        }
                    }
                }
            }
        }

        #pragma omp parallel for collapse(2) num_threads(32)
        for (int i=1;i<=nodesCount;++i){
            for (int j=1;j<=nodesCount;++j){
                if (diameter<FWDist[i][j]){
                    diameter=FWDist[i][j];
                }
            }
        }
    }
    //et = omp_get_wtime();
    //printf("run time: %lf s\n",et-st);
    printf("%d\n", diameter);
    return 0;
}

