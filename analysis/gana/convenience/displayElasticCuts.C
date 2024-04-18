#include <TCanvas.h>
#include <TPaveText.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include "../../../src/jsonmgr.C"

void parseAndDisplayCuts(const char* cuts) {
    // Convert the input string to a ROOT TString for easy manipulation
    TString strCuts(cuts);
    
    // Parse the conditions separated by "&&" into a vector of TStrings
    TObjArray* tokens = strCuts.Tokenize("&&");
    if (!tokens) return;
    
    // Create a canvas to display the conditions
    TCanvas* c1 = new TCanvas("c1", "Elastic Cuts Display", 800, 600);
    
    // Create a TPaveText to hold and display the conditions
    TPaveText* pt = new TPaveText(0.1, 0.1, 0.9, 0.9);
    pt->AddText("Elastic Cuts:");
    pt->AddText(" "); // Adding an empty line for better readability
    
    // Loop over the parsed conditions and add each to the TPaveText
    for (int i = 0; i < tokens->GetEntries(); ++i) {
        TObjString* token = (TObjString*)tokens->At(i);
        TString condition = token->GetString();
        //condition.ReplaceAll("abs(", ""); // Optional: Simplify display by removing 'abs('
        //condition.ReplaceAll(")", "");    // Optional: Simplify display by removing ')'
        pt->AddText(condition.Data());
    }
    
    // Draw the TPaveText on the canvas
    pt->Draw();
    
    // Update the canvas to display the changes
    c1->Update();

    // Clean up
    delete tokens; // Free memory used by the TObjArray
}

void displayElasticCuts(int kine = 9, int pass = 2, int mag = 70) {

  JSONManager *jmgr = new JSONManager("../../../config/syst.json");

  std::string cuts = jmgr->GetValueFromSubKey_str( Form("post_cuts_p%d",pass), Form("sbs%d_%d",kine,mag) );


  parseAndDisplayCuts(cuts.c_str());
}
