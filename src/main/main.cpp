#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "proc/cwdtw.h"
#include "util/exception.h"
#include "kmer/kmer_index.h"
#include "wavelib.h"
#include <malloc.h>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;


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
		bound = genomes.size()-4;//genomes.size()%5;
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
		bound = genomes.size()-5;//genomes.size()%5;
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
				signals.push_back(5.7*100+14);
			}
		}
	}

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
	opts.kmer    = 1;       //-> [1] to use 6mer; 0 to use 5mer
	

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
	
	long npyr = opts.level;          // default: 3
	double scale0 = opts.scale0;	  // default: sqrt(2)
	double dscale = 1;              // default: 1

	g::cwdtw::CWTAnalysis(reference, rcwt, scale0, dscale, npyr);	
	g::cwdtw::CWTAnalysis(peer, pcwt, scale0*alpha, dscale, npyr);

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
	g::cwdtw::MultiLevel_WaveletDTW(reference, peer, rcwt, pcwt, cosali, opts.radius, opts.test, opts.mode, &tdiff);
	if(opts.verbose ==1){
		EX_TRACE("Average Deviation (%.1lf/%ld=%.3lf)\n", tdiff, cosali.size(), tdiff/cosali.size());
	}

	//------ 5.1 generate final boundary -------//
	std::vector<std::pair<long,long> > bound;
	g::cwdtw::BoundGeneration(cosali, opts.neib, bound, opts.mode, 1);


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

