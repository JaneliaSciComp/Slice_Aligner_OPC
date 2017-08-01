# Slice_Aligner_OPC
This is a Fiji plugin for aligning the 2D + time sequences of bioimaging data.
This plugin aligns all slices within a 3D stack to a single template slice. The template slice is needed to be chosen from a stack. <br><br>
Written by Hideo Otsuna (HHMI Janelia, Scientific Computing Software, Senior Scientist)

<br>
<div align="right">
<a href="https://www.janelia.org/"><img src="https://github.com/JaneliaSciComp/Slice_Aligner_OPC/blob/master/maual/jrc_logo_180x40.png" alt="Link to Janelia"></a></div>
<br>
# **Manual**
<br>

<img src="https://github.com/JaneliaSciComp/Slice_Aligner_OPC/blob/master/maual/GUI.jpg" width="500" height="355" alt="GUI">
<br>
<br>

**Template slice:**
  * Current slice: Use current slice as a template.
  * First slice: Slice number 1 for the template.
  * input slice number: Typing window will appear after this GUI, then type slice number of the template in the window.

<br>

**+ rotation angle:**
  * Maximum positive rotation angle to search image overlap. (Less number is faster)
 <br>
 
**- rotation angle:**
  * Maximum negative rotation angle to search image overlap. (Less number is faster)
<br>

**Overlap percentage:**
  * Template and sample overlap percentage. 100 means 0 pixel shifting to search overlap, 50 means +-50% image shifting. (Higher number is faster) 
<br>

**Parallel Threads:**
  * Number of CPU for this plugin usage.
<br>

**Rotation increment:**
  * Increment angle for the rotation searching.
<br>

**Reference channel:**
  * C1, C2, All channels (C1 + C2).
<br>

**70% background subtraction:**
  * Subtract background before the template/sample matching calculation. This function will help the accuracy of the alignment with high-background sample.

<br>

