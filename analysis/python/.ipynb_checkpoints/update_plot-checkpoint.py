import ROOT
from ipywidgets import interact, FloatSlider
import numpy as np

# Load your ROOT file
file = ROOT.TFile.Open("/lustre19/expphy/volatile/halla/sbs/seeds/parse/parse_mc_sbs4_30p_barebones_alt.root")
tree = file.Get("P")

# Define a function to draw the plot based on the sliders' positions
def update_plot(x_cut, y_cut):
    # Clear the current canvas
    ROOT.gROOT.FindObject('c1').Clear()

    # Define the cut string based on slider values
    cut_string = f"mc_weight_norm*(bb_ps_e>{x_cut} && W2<{y_cut})"
    
    # Use TTree::Draw with the cut string. Adjust variable names and binning as needed.
    tree.Draw(f"dx:W2>>hist(100,-10,10,100,0,10)", cut_string, "colz")
    
    # Draw the updated histogram
    ROOT.gPad.Update()

# Create interactive sliders
interact(update_plot, 
         x_cut=FloatSlider(value=1, min=0, max=10, step=0.1, description='X Cut:'),
         y_cut=FloatSlider(value=1, min=0, max=10, step=0.1, description='Y Cut:'))

# Initial plot
c1 = ROOT.TCanvas('c1', '', 800, 600)
update_plot(1, 1)
