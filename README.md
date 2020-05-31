# <font color="CC0000">Claustrum Classifier</font>

Automated electrophysiology based classification of claustrum neurons

### Aditya Nair, Martin Graf, George Augustine<sup><a href="#fn1" id="ref1">1</a></sup>

The claustrum is a mysterious nucleus of the brain that regulates diverse behaviour from sleep to attention and possibly even consciousness! At the heart of such myriad behaviours might be the presence of heterogenous cell-types in the claustrum. 

The Augustine lab has recently performed an extensive classification of claustrum neurons based on intrinsic electrophysiological properties and identified at least 5 such cell-types<sup><a href="#fn2" id="ref2">2</a></sup>! The Claustrum Classifier allows anyone to use this scheme by providing software tools for automated extraction of cellular properties from ex-vivo electrophysiological data and classification using a trained neural net that can distinguish the cell-types identified in our publication.

## About Claustrum Classifier

> Claustrum Classifier is written in the R programming language and uses the Shiny framework for its GUI. It uses several packages internally including abf2, ggplot2, peakdet  

## Using Claustrum Classifier

# Use our online webapp!

> The easiest way to use the software is at the following link as a webapp: [Claustrum Classifier](https://claustrum.shinyapps.io/online/). A step-by-step walkthrough is given our Wiki page on uploading data and interpreting results.

# Install and use locally 

> If you wish to use the software on your computer, you would first need to use the script Install_Packages_CLA.R to install the various packages used. After this, simply perform 'runApp('CLA_Classifier.R')'

## Ease of Use 

>SilberReconstructor acts as a pipeline that allows for the complete automated reconstruction of neuronal morphology. It presents a simple to use interface that allows users to proceed with analysis without any domain-specific knowledge in image processing.

## License

Claustrum Classifier is licensed through Nanyang Technological University; redistribution and use for academic and other non-commercial purposes, with or without modification, are permitted provided that conditions of the license are met

--- 
 

### For any queries, please contact [Aditya Nair](adi.nair@caltech.edu)

 

<sup id="fn1">1. [Augustine Laboratory, Lee Kong Chian School of Medicine, Nanyang Technological University](http://www.lkcmedicine.ntu.edu.sg/aboutus/Faculty-and-Staff/Pages/George-Augustine.aspx)<a href="#ref1" title="Jump back to footnote 1 in the text.">↩</a></sup> 

<sup id="fn1">2. [Graf, Nair, Wong, Tang and Augustine, Identification of mouse claustral neuron types based on their intrinsic electrical properties, Submitted for publication](https://www.abstractsonline.com/pp8/#!/4376/presentation/33214)<a href="#ref2" title="Jump back to footnote 2 in the text.">↩</a></sup> 
