# pLM4ACE

## The implementation of the paper pLM4ACE: A protein language model-based deep learning predictor for screening peptide with high antihypertensive activity
Web server site: [https://sqzujiduce.us-east-1.awsapprunner.com]

The original scripts for webserver development are available at [https://github.com/dzjxzyd/LM4ACE_webserver]

#### Notice: pLM4ACE ONLY freely available for academic research; for commercial usage, please contact us

If the contents are useful to you, Please kindly Star it and Cite it. Please cite: pLM4ACE: A protein language model-based deep learning predictor for screening peptide with high antihypertensive activity



## Requirements
The majoy dependencies used in this project are as following:
```
Flask==2.1.0  # for web server 
pandas==1.3.5
numpy==1.21.6
fair-esm==2.0.0
scikit-learn===1.0.2
torch==1.12.0+cpu
pickle==4.0  # for web server 
gunicorn # for web server 
```
More detailed python libraries used in this project are referred to requirements.txt. All the implementation can be down in Google Colab and all you need is just a browser and a google account. Install all the above packages by !pip install package_name==2.0.0

## Usage
Notice: all my dataset use 0 and 1 to represent positive and negative, respectively. Again, 0 is positive and 1 is negative.

## embedding approach references
ESM-2 https://github.com/facebookresearch/esm
iFeatureOmega https://github.com/Superzchen/iFeatureOmega-CLI
