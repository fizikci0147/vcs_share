//read in cmr
void ImportData(vector<double> &Q2,vector<double> &GE,vector<double> &E1, vector<double> &E2, vector<string> &author,const char*filename )
{
    double iQ2,iGE,iE1,iE2;
    string iauthor;
    ifstream infile(filename);
    if(infile.fail()){
        cout<<filename <<" input file doesn't exist!"<<endl;
        exit(1);
    }else{
        while(!infile.eof()){
            infile.ignore(272,'\n');
            infile>>iQ2>>iGE>>iE1 >>iE2>>iauthor;
            Q2.push_back(iQ2);
            GE.push_back(iGE);
            E1.push_back(iE1);
            E2.push_back(iE2);
            author.push_back(iauthor);
        }
        infile.close();
        Q2.pop_back();
        GE.pop_back();
        E1.pop_back();
        E2.pop_back();
        author.pop_back();
    }
}
