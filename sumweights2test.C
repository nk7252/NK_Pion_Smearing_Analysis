

void sumweights2test() {

    TH1 *h1= new TH1D("h1","unweighted version test", 64, 0, 64);
    TH2 *h2= new TH2D("h2","2D gauss, unweighted version test", 64, 0, 64, 64, 0, 64);
    TH2 *h3= new TH2D("h3","2D gauss, weighted version test", 64, 0, 64, 64, 0, 64);
    auto hprof = new TProfile("hprof","Profile test",64, 0, 64, 0, 64);
    auto hprof2 = new TProfile("hprof2","Profile test2",64, 0, 64, 0, 64);
    //h3->Sumw2();
    hprof2->Sumw2();
    //TF1 *f1 = new TF1();
    TRandom3 gen;
    double x, y, w;
    for (int i=0;i<1000000; i++){
        x=gen.Gaus(30.,5.);
        y=gen.Gaus(30.,4.);
        w=pow(sqrt(pow(x,2)+pow(y,2)),5);
        h1->Fill(x);
        h2->Fill(x,y);
        h3->Fill(x,y,w);
        hprof->Fill(x,y,w);
        hprof2->Fill(x,y,w);
    }
    //h1->Draw("COLZ");
    //h2->Draw("COLZ");
    auto c1= new TCanvas("c1","test canvas",200,10,1200,300);
    c1->Divide(3,1);
    c1->cd(1);
    h3->Draw("COLZ");
    c1->cd(2);
    hprof->Draw("COLZ");
    c1->cd(3);
    hprof2->Draw("COLZ");
}