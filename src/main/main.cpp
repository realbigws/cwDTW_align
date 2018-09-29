#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "proc/proc.h"
#include "util/exception.h"
#include "kmer/kmer_index.h"
#include "wavelib.h"
#include <malloc.h>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace g::proc;

//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}
//---------------------- utility ------//over



//--------- pore_models: from genome to expected signal ----------//
void Genomes2SignalSequence(const std::vector<char> &genomes, 
	std::vector<int>& index, std::vector<double>& signals, 
	int scale, int FIVE_or_SIX, int ZSCO_or_NOT)
{
	long bound;
	if(FIVE_or_SIX==0) //-> 5mer model
	{
		g::Mer2Signal::Genome2Index_5mer(genomes, index);
		bound = genomes.size()-5;//genomes.size()%5;
		signals.assign(bound*scale,0);
		long cur=0;
		for(long i = 0; i < bound; i++){
			double sigval;
			if(index[i]<0)sigval=0;
			else{
				sigval = g::Mer2Signal::AvgSignalAt_5mer(index[i]);
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					sigval = (sigval-90.208351)/12.832660;
				}
				else               //-> use original int value
				{
					sigval = (int)(5.7*sigval+14);
				}
			}
			for(int c = scale; c--;){
				signals[cur]=sigval;
				cur++;
			}
		}
	}
	else               //-> 6mer model
	{
		g::Mer2Signal::Genome2Index_6mer(genomes, index);
		bound = genomes.size()-6;//genomes.size()%5;
		signals.assign(bound*scale,0);
		long cur=0;
		for(long i = 0; i < bound; i++){
			double sigval;
			if(index[i]<0)sigval=0;
			else{
				sigval = g::Mer2Signal::AvgSignalAt_6mer(index[i]);
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					sigval = (sigval-90.208199)/12.868652;
				}
				else               //-> use original int value
				{
					sigval = (int)(5.7*sigval+14);
				}
			}
			for(int c = scale; c--;){
				signals[cur]=sigval;
				cur++;
			}
		}
	}

	//---- tail five_mer ------//
	for(long i = bound; i < genomes.size(); i++){
		for(int c = scale; c--;){
			if(ZSCO_or_NOT==1) //-> transfer to Zsco
			{
				signals.push_back(0);
			}
			else
			{
				signals.push_back(100);
			}
		}
	}

}

//--------------- continuous wavelet transform (CWT) analysis -----------------//
/** @scale0: level0 pyramind scale;  @dscale: scale_i = scale0*(2^{i*dsacle} ); @npyr: total number of pyramind*/
void CWTAnalysis(const std::vector<double>& raw, std::vector<std::vector<double> >& output, 
	double scale0, double dscale, int npyr)
{
	const double* sigs = &raw[0];		//sst_nino3.dat
	cwt_object wt;

	size_t N = raw.size();
	double dt = 1;//2;		//sample rate	>  maybe we should use 2?

	wt = cwt_init("dog", 2.0, N, dt, npyr);	//"morlet", "dog", "paul"
	setCWTScales(wt, scale0, dscale, "pow", 2.0);
	cwt(wt, sigs);

	output.resize(npyr);
	for(size_t k = npyr; k--;){
		int idx = npyr-k-1;
		
		output[idx].resize(raw.size());
		size_t offset = k*raw.size();
		for(size_t i = 0; i < output[idx].size(); i++){
			output[idx][i] = wt->output[i+offset].re;
		}
	}

	cwt_free(wt);
}

//----------------- boundary generation for constrained dynamic time warping (cDTW) -------------//
void BoundGeneration(std::vector<std::pair<long,long> >& cosali, 
	long neib, std::vector<std::pair<long,long> >& bound, int mode,
	int RENMIN_or_SHENG=0)
{

//if(mode!=-1) //-> mode = -1 means Renmin mode
vector<pair<long,long> > cosali_=cosali;
if(mode==1) //-> use partial-diagonol alignment
{
	//-------- generate partial-diagonol alignment -------//
	vector<pair<long,long> > cosali2;
	cosali2.push_back(cosali[0]);
	for(long i = 1; i < cosali.size(); i++){
		if(cosali[i].first != cosali2[cosali2.size()-1].first){
			cosali2.push_back(cosali[i]);
		}
		else{
			cosali2[cosali2.size()-1].second = cosali[i].second;
		}
	}

	cosali_.clear();
	for(long i = 1; i < cosali2.size(); i++)
	{
		long pre_idx = cosali2[i-1].first, cur_idx = cosali2[i].first;
		long pre_anchor = cosali2[i-1].second, cur_anchor = cosali2[i].second;
		double anchor_diff = 1.0*(cur_anchor-pre_anchor)/(double)(cur_idx-pre_idx);
		for(long k = pre_idx, count = 0; k < cur_idx; k++, count++)
		{
			long mid = pre_anchor+(long)(count*anchor_diff);  //assign point relationship
			if(mid<pre_anchor)mid=pre_anchor;
			if(mid>cur_anchor)mid=cur_anchor;
			cosali_.push_back(make_pair(k,mid));
		}
	}
	cosali_.push_back(cosali2[cosali2.size()-1]);
}

	//----> output pre-bound alignnment -------
//	{
//		for(int i = 0; i < cosali_.size(); i++)
//		{
//			printf("%d -> %d  %d \n",i,cosali_[i].first,cosali_[i].second);
//		}
//	}


	//---------- use block bound ------------//start
	//-> get signal length
	long moln1=cosali_[cosali_.size()-1].first+1;
	long moln2=cosali_[cosali_.size()-1].second+1;
	//-> renmin align to sheng style
	std::vector<std::pair<long,long> > align_sheng;
	Renmin_To_Sheng_align(moln1,moln2,cosali_,align_sheng);
	//-> get bound in sheng style
	std::vector<std::pair<long,long> > bound_sheng;
	From_Align_Get_Bound(moln1,moln2,align_sheng,bound_sheng,neib);
	//-> transfer bound to renmin style
	Sheng_To_Renmin_bound(moln1,moln2,bound_sheng,bound);
	//----- maybe useful -----//
	bound[0].first = 0;
	bound[bound.size()-1].first = bound[bound.size()-2].first;
	bound[bound.size()-1].second = cosali[cosali.size()-1].second;


	if(RENMIN_or_SHENG==1) //-> use Sheng's bound definition
	{
		bound=bound_sheng;
	}

	//----> output post-bound alignnment -------
//	{
//		for(int i = 0; i < bound.size(); i++)
//		{
//			printf("%d -> %d  %d \n",i,bound[i].first,bound[i].second);
//		}
//		exit(-1);
//	}

	return;
	//---------- use block bound ------------//over


/*
	int radius=neib1;
	std::vector<std::pair<int,int> > cosali2;
	
	cosali2.push_back(cosali[0]);
	for(int i = 1; i < cosali.size(); i++){
		if(cosali[i].first != cosali2[cosali2.size()-1].first){
			cosali2.push_back(cosali[i]);
		}
		else{
			cosali2[cosali2.size()-1].second = cosali[i].second;
		}
	}
	
	cosali = cosali2;
	
	bound.resize(cosali[cosali.size()-1].first+1);
	bound.assign(cosali[cosali.size()-1].first+1, std::make_pair(-1,-1));
	
	for(int i = 1; i < cosali.size(); i++){
		int pre_idx = cosali[i-1].first, cur_idx = cosali[i].first;
		int pre_anchor = cosali[i-1].second, cur_anchor = cosali[i].second;
		
		double anchor_diff;
		
		anchor_diff = (cur_anchor-pre_anchor)/(cur_idx-pre_idx);
		
		int neighbor = int(anchor_diff+0.5)+radius;
		
		for(int i = pre_idx, count = 0; i < cur_idx; i++, count++){
			int mid = pre_anchor+std::round(count*anchor_diff);		//assign point relationship
			bound[i].first = mid-neighbor < 0 ? 0 : mid-neighbor;
			bound[i].second = mid+neighbor > cosali[cosali.size()-1].second ? cosali[cosali.size()-1].second : mid+neighbor;
		}
	}
	
// 	for(int i = 0; i < bound.size(); i++){
// 		std::cout<<bound[i].first<<" "<<bound[i].second<<std::endl;
// 	}
	
	bound[0].first = 0;
	bound[bound.size()-1].first = bound[bound.size()-2].first;
	bound[bound.size()-1].second = cosali[cosali.size()-1].second;
*/
}

//====================== continous wavelet dynamic time warping (cwDTW) ========================//
void MultiLevel_WaveletDTW(std::vector<double>& in1, std::vector<double>& in2,
	std::vector<std::vector<double> >& sig1, std::vector<std::vector<double> >& sig2, 
	std::vector<std::pair<long,long> >& alignment, long radius, int test, int mode, 
	double* totaldiff = 0)
{
	double tdiff;
	std::vector<std::pair<long, double> > sig1peaks, sig2peaks;
	double length1 = sig1[0].size(), length2 = sig2[0].size();

	long tot_size=sig1.size();
	for(long k = 0; k < tot_size; k++)
	{
		//------ peakpick CWT signal -------------//
		g::proc::PeakPick(sig1[k], sig1peaks);
		g::proc::PeakPick(sig2[k], sig2peaks);
//		std::cout<<sig1peaks.size()<<"\t"<<sig2peaks.size()<<std::endl;
		std::vector<double> peak1(sig1peaks.size());
		std::vector<double> peak2(sig2peaks.size());


		//-------- use peak picking result ---------//
		for(long i = sig1peaks.size(); i--;){
			peak1[i] = sig1peaks[i].second;
		}
		for(long i = sig2peaks.size(); i--;){
			peak2[i] = sig2peaks[i].second;
		}

		//----- apply DTW or cDTW dependent on k-th level -------//
		if(k == 0){
			tdiff = g::proc::DynamicTimeWarping(peak1, peak2, alignment);
		}
		else{
			//----- ReMapIndex_partI (map ground level upon k-th level) -----//
			long c = 0;
			for(long i = 0; i < alignment.size(); i++){
				while(sig1peaks[c].first < alignment[i].first){
					c++;
				}
				alignment[i].first = c;	
			}
			
			long d = 0;
			for(long i = 0; i < alignment.size(); i++){
				while(sig2peaks[d].first < alignment[i].second){
					d++;
				}
				alignment[i].second = d;	
			}
			//----- cDWT (constrained DWT) -------//
			std::vector<std::pair<long,long> > bound;
			long neib=radius*pow(2,tot_size-k); //adaptive radius
			BoundGeneration(alignment, neib, bound, mode);
			tdiff = g::proc::BoundDynamicTimeWarping(peak1, peak2, bound, alignment);
		}
		//----- ReMapIndex_partII (map k-th level back to ground level) -----//
		for(long i = alignment.size(); i--;){
			alignment[i].first = sig1peaks[alignment[i].first].first;
			alignment[i].second = sig2peaks[alignment[i].second].first;
		}

	}
	
	if(totaldiff){
		*totaldiff = tdiff;
	}
}

//------------------ write alignment to file ----------------------//
void WriteSequenceAlignment(const char* output,
        const std::vector<double>& reference, const std::vector<double>& peer,
        vector<pair<long,long> >& alignment, int swap)
{
	vector <std::string> tmp_rec;
	double diff;
	for(long i = 0; i < alignment.size(); i++)
	{
		long pos1=alignment[i].second;
		long pos2=alignment[i].first;
		//----- output to string ----//
		std::ostringstream o;
		if(swap==1)
		{
			o<<setw(10)<<alignment[i].second<<" "<<setw(9)<<alignment[i].first<<" | ";
			o<<setw(15)<<peer[alignment[i].second]<<", "<<setw(15)<<reference[alignment[i].first];
		}
		else
		{
			o<<setw(10)<<alignment[i].first<<" "<<setw(9)<<alignment[i].second<<" | ";
			o<<setw(15)<<peer[alignment[i].first]<<", "<<setw(15)<<reference[alignment[i].second];
		}
		o<<"          diff:"<<setw(15)<<(diff = std::fabs(reference[alignment[i].first]-peer[alignment[i].second]));
		//----- record string -----//
		std::string s=o.str();
		tmp_rec.push_back(s);
	}
	//----- output to file ------//
	FILE *fp=fopen(output,"wb");
	for(long i=0;i<tmp_rec.size();i++)fprintf(fp,"%s\n",tmp_rec[i].c_str());
	fclose(fp);
}


//================================== DNA_DynaProg_NanoPore_assess_FixFast =============================//

//=================== NanoPore DynaProg ===============//
//-> ATCG to digit
//-> see below link for base pairing:
//       http://www.biology-pages.info/B/BasePairing.html
//-> see below link for molecular weight:
//       http://www.genomics.agilent.com/files/Mobio/Nucleic%20Acids_Sizes_and_Molecular_Weights_2pgs.pdf
int DNA_To_Int(char c)
{
	switch(c)
	{
		case 'A': return -2;
		case 'T': return -1;
		case 'C': return  1;
		case 'G': return  2;
		return 0;
	}
}

//-------- DNA calc ------//
int DNA_Calc(char a,char b)
{
	int mat=1;
	int mis=-2;
	int s1=DNA_To_Int(a);
	int s2=DNA_To_Int(b);
	if(s1==0 || s2==0)return mis;
	if(a==b)return mat;
	else return mis;
}

//---- part 2. generate alignment from a given bound ------//
int Advance_Align_Dyna_Prog_Double_Bound(long n1,long n2,const vector <vector<double> > &score,
	double GAP_OPEN1,double GAP_EXT1,double GAP_OPEN2,double GAP_EXT2,
	double GAP_HEAD1,double GAP_TAIL1,double GAP_HEAD2,double GAP_TAIL2,
	vector<pair<long,long> > &bound, vector<pair<long,long> > & alignment,double &ali_sco)
{
	long i,j,k;
	//input
	long m = n1 + 1;  // +1 to account for the extra row,col in
	long n = n2 + 1;  // the DP matrices corresponding to gaps
	//const value
	const int _H_  = 0;
	const int _S_  = 1;
	const int _V_  = 2;

	//create D and M
	vector <vector <vector <int> > > D;      // the path (directions) matrix
	vector <vector <vector <double> > > M;   // the current scores (values) matrix
	D.resize(3);
	M.resize(3);
	//resize(m,n)
	long wsstart,wsend;
	long prev_start,prev_end;
	for (k = 0; k < 3; k++) 
	{
		D[k].resize(m);
		M[k].resize(m);
		//-> init
		for(i=0; i<m; i++)
		{
			//current bound
			wsstart=bound[i].first;
			wsend=bound[i].second;
			long bound_len=wsend-wsstart+1;
			D[k][i].resize(bound_len,-1);
			M[k][i].resize(bound_len,-1);
		}
	}
	//init()
	double WS_MIN=-1000000;
	D[_S_][0][0] = -1;
	D[_H_][0][0] = -1;
	D[_V_][0][0] = -1;
	M[_S_][0][0] = 0;
	M[_H_][0][0] = WS_MIN;
	M[_V_][0][0] = WS_MIN;
	//first line
	i=0;
	wsstart=bound[i].first;
	wsend=bound[i].second;
	for (j = 1; j <= wsend; j++) 
	{
		D[_S_][i][j-wsstart] = _H_;
		D[_H_][i][j-wsstart] = _H_;
		D[_V_][i][j-wsstart] = _H_;
		M[_S_][i][j-wsstart] = WS_MIN;
		M[_H_][i][j-wsstart] = j*GAP_HEAD2; //-(Params::GAP_OPEN + (j-1)*Params::GAP_EXT);
		M[_V_][i][j-wsstart] = WS_MIN;
	}
	prev_start=wsstart;
	prev_end=wsend;
	//fill(firstSeq, secondSeq, distFunc);
	double gap_open;
	double gap_ext;
	double v1,v2,v3;
	double dist;
	for (i = 1; i < m; i++) 
	{
		//current bound
		wsstart=bound[i].first;
		wsend=bound[i].second;
		//process
		for (j = wsstart; j <= wsend; j++) 
		{
			//------ 1. condition upper -----//
			if(j==n-1)
			{
				gap_open=GAP_TAIL1;
				gap_ext=GAP_TAIL1;
			}
			else
			{
				gap_open=GAP_OPEN1;
				gap_ext=GAP_EXT1;
			}
			//-- boundary --//(start)
			if(j>prev_end)
			{
				M[_V_][i][j-wsstart] = WS_MIN;
				D[_V_][i][j-wsstart] = -1;
			}
			else
			{
				if(j==0)
				{
					gap_ext=GAP_HEAD2;
					gap_open=GAP_HEAD2;
				}
				long curwsstart=bound[i-1].first;
				v1 = M[_V_][(i-1)][j-curwsstart] + gap_ext;
				v2 = M[_S_][(i-1)][j-curwsstart] + gap_open;
				v3 = M[_H_][(i-1)][j-curwsstart] + gap_open;
				M[_V_][i][j-wsstart] = std::max(v1, std::max(v2, v3));
				if (M[_V_][i][j-wsstart] == v1) D[_V_][i][j-wsstart] = _V_;
				else if(M[_V_][i][j-wsstart] == v2) D[_V_][i][j-wsstart] = _S_;
				else D[_V_][i][j-wsstart] = _H_;
			}
			//-- boundary --//(end)

			//------ 2. condition left -----//
			if(i==m-1)
			{
				gap_open=GAP_TAIL2;
				gap_ext=GAP_TAIL2;
			}
			else
			{
				gap_open=GAP_OPEN2;
				gap_ext=GAP_EXT2;
			}
			//-- boundary --//(start)
			if(j==wsstart)
			{
				M[_H_][i][j-wsstart] = WS_MIN;
				D[_H_][i][j-wsstart] = -1;
			}
			else
			{
				v1 = M[_H_][i][j-1-wsstart] + gap_ext;
				v2 = M[_S_][i][j-1-wsstart] + gap_open;
				v3 = M[_V_][i][j-1-wsstart] + gap_open;
				M[_H_][i][j-wsstart] = std::max(v1, std::max(v2, v3));
				if (M[_H_][i][j-wsstart] == v1) D[_H_][i][j-wsstart] = _H_;
				else if(M[_H_][i][j-wsstart] == v2) D[_H_][i][j-wsstart] = _S_;
				else D[_H_][i][j-wsstart] = _V_;
			}
			//-- boundary --//(end)

			//------ 3. condition diag -----//
			//-- boundary --//(start)
			if(j==prev_start || j>prev_end+1)
			{
				M[_S_][i][j-wsstart] = WS_MIN;
				D[_S_][i][j-wsstart] = -1;
			}
			else
			{
				long curwsstart=bound[i-1].first;
				long start=wsstart-1<0?0:wsstart-1;
				dist = score[i-1][j-1-start];
				v1 = M[_V_][(i-1)][j-1-curwsstart] + dist;
				v2 = M[_H_][(i-1)][j-1-curwsstart] + dist;
				v3 = M[_S_][(i-1)][j-1-curwsstart] + dist;
				M[_S_][i][j-wsstart] = std::max(v1, std::max(v2, v3));
				if (M[_S_][i][j-wsstart] == v3) D[_S_][i][j-wsstart] = _S_;
				else if (M[_S_][i][j-wsstart] == v1) D[_S_][i][j-wsstart] = _V_;
				else D[_S_][i][j-wsstart] = _H_;
			}
			//-- boundary --//(end)
		}
		//next bound
		prev_start=wsstart;
		prev_end=wsend;
	}

	//build(ali, firstSeq, secondSeq, distFunc);
	i = m-1;
	wsstart=bound[i].first;
	j = n-1;
	v1=M[_V_][i][j-wsstart];
	v2=M[_H_][i][j-wsstart];
	v3=M[_S_][i][j-wsstart];
	double maximal = std::max(v1, std::max(v2, v3));
	k = -1;
	if(v3==maximal)k = _S_;
	else if(v2==maximal)k = _H_;
	else k = _V_;
	//trace_back
	alignment.clear();
	long count = 0;
	long matches = 0;
	long cur_case=k;
	long pre_case;
	for(;;)
	{
		if(i==0||j==0)break;
		wsstart=bound[i].first;
		wsend=bound[i].second;
		pre_case=D[cur_case][i][j-wsstart];
		switch (cur_case)
		{
			case _S_:
				alignment.push_back(pair<long,long>(i,j)); 
				i--;
				j--;
				++matches;
				break;
			case _V_:
				alignment.push_back(pair<long,long>(i,-j)); 
				i--;
				break;
			case _H_:
				alignment.push_back(pair<long,long>(-i,j)); 
				j--;
				break;
			default:
				cout << "ERROR!! -> advance_global: invalid direction D[" << k << "](" << i << ", " << j << ") = " 
				<< D[k][i][j-wsstart] << endl;
				exit(-1);
		}
		cur_case=pre_case;
		count++;
	}
	while (j> 0) alignment.push_back(pair<long,long>(-i,j)),j--;
	while (i> 0) alignment.push_back(pair<long,long>(i,0)), i--;
	reverse(alignment.begin(), alignment.end());
	ali_sco=maximal;
	return matches;
}

//============= Output Alignment ========//
//-> fasta_output_simp
void FASTA_Output_Simp(FILE *fp,string &nam1,string &nam2,
	const char *ami1,const char *ami2,
	vector<pair<long,long> > &alignment,
	string &assess)
{
	//alignment->output
	//output
	char c;
	long i;
	long ii,jj;
	long size=alignment.size();
	//--> output seq1
	fprintf(fp,">%s\n",nam1.c_str());
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		if(ii<=0)fprintf(fp,"-");
		else
		{
			c=ami1[ii-1];
			fprintf(fp,"%c",c);
		}
	}
	fprintf(fp,"\n");
	//--> output assess
	if(assess!="")
	{
		fprintf(fp,"%s\n",assess.c_str());
	}
	//--> output seq2
	fprintf(fp,">%s\n",nam2.c_str());
	for(i=0;i<size;i++)
	{
		jj=alignment[i].second;
		if(jj<=0)fprintf(fp,"-");
		else
		{
			c=ami2[jj-1];
			fprintf(fp,"%c",c);
		}
	}
	fprintf(fp,"\n");
}



//----------------- main -------------------//
int main(int argc, char **argv)
{
	struct options opts;
	opts.radius  = 15;
	opts.level   = 4;
	opts.scale0  = sqrt(2);
	opts.neib    = 25;
	opts.verbose = 0;       //-> [0] no verbose; 1 verbose
	opts.test    = 0;       //-> [0] not use test mode; 1 equal_ave, 2 peak_ave, 3 Fast_DTW
	opts.mode    = 0;       //-> [0] block bound; 1 diagonol bound
	opts.kmer    = 0;       //-> [0] to use 5mer; 1 to use 6mer
	

	//----- parse arguments -----//
	if(GetOpts(argc, argv, &opts) < 0){
		EX_TRACE("**WRONG INPUT!**\n");
		return -1;
	}

	std::string input1=opts.input;
	std::string input2=opts.peer;
	std::string output=opts.output;
	if(input1=="" || input2=="")
	{
		fprintf(stderr,"input1 or input2 is NULL \n");
		return -1;
	}


	//======================= START Procedure ===================================//

	//=========================================//
	//----- 1. read genome translated signal -----//
	std::vector<char> genomes1;   //genome sequence
	if(!g::io::ReadATCG(opts.input, genomes1)){
		EX_TRACE("Cannot open %s.\n", opts.input);
		return -1;
	}
	//----- 1.1 pore_model: transform genome sequence to expected signal -------//
	std::vector<int> refer_orig;
	std::vector<double> reference;  //reference: genome signal
	Genomes2SignalSequence(genomes1, refer_orig, reference, 1, opts.kmer, 1);
	//----- 1.2 get input name1 --------//
	std::string genom_name_orig=opts.input;
	std::string genom_name;
	getBaseName(genom_name_orig,genom_name,'/','.');


	//========================================//
	//----- 2. read genome translated signal -----//
	std::vector<char> genomes2;   //genome sequence
	if(!g::io::ReadATCG(opts.peer, genomes2)){
		EX_TRACE("Cannot open %s.\n", opts.peer);
		return -1;
	}
	//----- 2.1 pore_model: transform genome sequence to expected signal -------//
	std::vector<int> peer_orig;
	std::vector<double> peer;  //reference: genome signal
	Genomes2SignalSequence(genomes2, peer_orig, peer, 1, opts.kmer, 1);
	//----- 2.2 get input name1 --------//
	std::string signal_name_orig=opts.peer;
	std::string signal_name;
	getBaseName(signal_name_orig,signal_name,'/','.');


//exit(-1);    //-> 2.038047 seconds

	//==================================================//
	//------3. process initial input signals ----------//

	//----- 3.1 Zscore normaliza on both signals -----//
//	g::proc::ZScoreNormalize(reference);
//	g::proc::ZScoreNormalize(peer);

	//----- 3.2 calculate length ratio between input signals -----//	
//	double alpha = (double)peer.size()/reference.size();
	double alpha = 1;


	//====================================================//
	//----- 4. continous wavelet transform --------------//
	std::vector<std::vector<double> > rcwt, pcwt;

	if(opts.verbose ==1){
		EX_TRACE("CWT Analysis...\n");
	}
	
	int npyr = opts.level;          // default: 3
	double scale0 = opts.scale0;	  // default: sqrt(2)
	double dscale = 1;              // default: 1

	CWTAnalysis(reference, rcwt, scale0, dscale, npyr);	
	CWTAnalysis(peer, pcwt, scale0*alpha, dscale, npyr);

	//------ 4.1 Zscore normaliza on both CWT signals -----//	
	//if multiscale is used, pyr logical should be added.
//	for(int i = 0; i < npyr; i++){
//		g::proc::ZScoreNormalize(rcwt[i]);
//		g::proc::ZScoreNormalize(pcwt[i]);
//	}

//exit(-1);   //-> 2.162879 seconds

	//============================================//
	//------ 5. multi-level WaveletDTW ----------//
	std::vector<std::pair<long,long> > cosali;
	double tdiff;

	if(opts.verbose ==1){
		EX_TRACE("Coarse Alignment...\n");
	}
	MultiLevel_WaveletDTW(reference, peer, rcwt, pcwt, cosali, opts.radius, opts.test, opts.mode, &tdiff);
	if(opts.verbose ==1){
		EX_TRACE("Average Deviation (%.1lf/%ld=%.3lf)\n", tdiff, cosali.size(), tdiff/cosali.size());
	}

	//------ 5.1 generate final boundary -------//
	std::vector<std::pair<long,long> > bound;
	BoundGeneration(cosali, opts.neib, bound, opts.mode, 1);


//exit(-1);   //-> 3.386344 seconds


	//------ 5.2 generate final alignment via cDTW ------//
//	std::vector<std::pair<int,int> > alignment;
//	tdiff = g::proc::BoundDynamicTimeWarping(reference, peer, bound, alignment);  //-> restrict version !!!
//	fprintf(stderr,"%s %s %lf %d %lf\n",signal_name.c_str(),genom_name.c_str(),tdiff,alignment.size(),tdiff/alignment.size());


	//=================================================//
	//------ 6. output final alignment to file -------//
//	if(output!="")
//		WriteSequenceAlignment_nano(opts.output, reference, peer, refer_orig, peer_orig, alignment, swap, genome_str);





	//======================== DNA_DynaProg_NanoPore_assess_FixFast ======================//
	{
		//-> 1. read fasta sequence
		std::string seqres1(genomes1.begin(), genomes1.end() );
		std::string seqres2(genomes2.begin(), genomes2.end() );
		long n1=seqres1.length();
		long n2=seqres2.length();

		//-> 2. dynamic programming alignment
		vector <vector <double> > BOUND_score;
		BOUND_score.resize(n1);
		for(long i=0;i<n1;i++)
		{
			//current bound
			long wsstart=bound[i+1].first-1;
			long wsend=bound[i+1].second-1;
			//process
			long start=wsstart<0?0:wsstart;
			long end=wsend>=n2?n2-1:wsend;
			for (long j = start; j <= end; j++) 
			{
				double score=DNA_Calc(seqres1[i],seqres2[j]);
				BOUND_score[i].push_back(score);
			}
		}
		double sco;
		vector<pair<long,long> > WWW_alignment;
		long matchs=Advance_Align_Dyna_Prog_Double_Bound(n1,n2,BOUND_score,-2,-1,-2,-1,0,0,0,0,
			bound,WWW_alignment,sco);

		//-- remove head and tail gap --//
		long iden=0;
		long igap=0;
		long ilen=0;
		long ogap=0;
		long mism=0;
		long first=1;
		string assess="";
		{
			//--> determine start
			long start=-1;
			for(long i=0;i<WWW_alignment.size();i++)
			{
				long ii=WWW_alignment[i].first;
				long jj=WWW_alignment[i].second;
				if(ii>0 && jj>0)
				{
					start=i;
					break;
				}
			}
			//--> determine end
			long end=-1;
			for(long i=WWW_alignment.size()-1;i>=0;i--)
			{
				long ii=WWW_alignment[i].first;
				long jj=WWW_alignment[i].second;
				if(ii>0 && jj>0)
				{
					end=i;
					break;
				}
			}
			//--> final check
			if(start==-1 || end==-1)
			{
				fprintf(stderr,"BAD HERE !!! Null Position Aligned !! \n");
				exit(-1);
			}
			//--> determine
			for(long i=0;i<start;i++)assess.push_back('.');
			for(long i=start;i<=end;i++)
			{
				long ii=WWW_alignment[i].first;
				long jj=WWW_alignment[i].second;
				if(ii>0 && jj>0)
				{
					if(seqres1[ii-1]==seqres2[jj-1])
					{
						assess.push_back('|');
						iden++;
					}
					else
					{
						assess.push_back('X');
						mism++;
					}
					//-- count gap open --//
					first=1;
				}
				else
				{
					assess.push_back('-');
					igap++;
					//-- count gap open --//
					if(first==1)
					{
						first=0;
						ogap++;
					}
				}
				ilen++;
			}
			for(long i=end+1;i<WWW_alignment.size();i++)assess.push_back('.');
		}

		//-> 3. output alignment
		string nam1=genom_name;
		string nam2=signal_name;
		if(output!="")
		{
			string ali_out=output;
			FILE *fp=fopen(ali_out.c_str(),"wb");
			FASTA_Output_Simp(fp,nam1,nam2,seqres1.c_str(),seqres2.c_str(),WWW_alignment,assess);
			fclose(fp);
		}
		printf("%s %s -> %d -> %d %d %d -> %lf %lf -> %d/%d(%lf) | %d/%d(%lf) -> %d/%d(%lf) | %d/%d(%lf) \n",
			nam1.c_str(),nam2.c_str(),(long)sco,
			iden,seqres1.length(),seqres2.length(),
			1.0*iden/seqres1.length(),1.0*iden/seqres2.length(),
			iden,ilen,1.0*iden/ilen,mism,ilen,1.0*mism/ilen,
			igap,ilen,1.0*igap/ilen,ogap,ilen,1.0*ogap/ilen);
	}

	//----- exit -----//	
	return 0;
}

