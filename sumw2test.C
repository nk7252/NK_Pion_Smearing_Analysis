void sumw2test()  {

  TH1F* hi = new TH1F("i", "i", 10, 0, 10);
  hi->Fill(3);
  std::cout << "sqrt(n) error " << hi->GetBinError(4) << std::endl;

  TH1F* hn =  new TH1F("n", "n", 10, 0, 10);
  hn->Fill(3, .5);
  hn->Fill(3, .5);
  hn->Sumw2(kFALSE);
  std::cout << "without SumW2 " << hn->GetBinError(4) << " uh oh" << std::endl;

  TH1F* hw = new TH1F("w", "w", 10, 0, 10);
//  hw->Sumw2(kTRUE); // not needed as Sumw2 is on by default
  //hw->Sumw2(kTRUE);
  hw->Fill(3, .5);
  hw->Fill(3, .5);
  std::cout << "with SumW2 " << hw->GetBinError(4) << std::endl;
}