# DeTox: Predicting Developmental Toxicity Using In-Silico Methods

Welcome to DeTox, an in-silico tool designed to replace animal testing by leveraging computational methods to predict developmental toxicity. This repository includes code for data ingestion, preprocessing, model inference, and a user-friendly web interface.

## Prerequisites
	•	Python ≥ 3.8
	•	Git

## Installation

Clone the repository and install dependencies:

```
git clone https://github.com/your-username/dev-tox.git
cd dev-tox
pip install -r requirements.txt
```
## Usage

The scraping tool is implemented in the dev_tox/scraper.py file.

## Running the Scraper

The scraper retrieves compound information from PubMed for downstream model processing.
```
1. Open python Webscraper/Webscraper_Pubmed.ipynb
2. Add compounds of interest and keywords to the input .csv file
3. Receive results in the Scrapper_results folder.
```

## Running the Web App

The online web tool can be accessed here:

# https://detox.mml.unc.edu/

The Webserver/.. folder contains a Flask application that allows users to interact with DeTox via a web interface.

### Locally running the Webserver

cd Webserver
python app.py

Visit http://127.0.0.1:5000 in your browser.


## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Reference
```
DeTox: an In-Silico Alternative to Animal Testing for Predicting Developmental Toxicity Potential. Ricardo Scheufen Tieghi, Marielle Rath, José Teófilo Moreira-Filho, James Wellnitz, Holli-Joi Martin, Kathleen Gates, Helena T. Hogberg, Nicole Kleinstreuer, Alexander Tropsha, Eugene N. Muratov.

```