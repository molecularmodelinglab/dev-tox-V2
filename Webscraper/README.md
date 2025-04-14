
# 🔧 Prerequisites
```
Ensure the following Python libraries are installed:

pip install pandas beautifulsoup4 requests openpyxl
```


⸻

# 📁 Input File Format
```
Prepare an Excel file named keywords_example.xlsx with a column titled "Compound_name" containing the list of compounds you wish to search on PubMed.

Example:

Compound_name

aspirin
paracetamol
ibuprofen
```


⸻

# 🚀 How to Run
```
	1.	Place your keywords_example.xlsx in the same directory as the script.
	2.	Execute the script in your terminal or IDE:

Execute the file "Webscraper_Pubmed.ipynb

The script will scrape the first 3 search results per keyword from PubMed and stop after your list finished or up to 600 compounds.
```
⸻

### 🧠 What the Script Does
```
1. Scrapes PubMed:
	•	Searches each compound name on PubMed.
	•	Extracts:
        •	Title
        •	Authors
        •	Journal
        •	Abstract

2. Curates Articles:

Filters and retains articles meeting either of the following criteria:
	•	Human-related: Mentions terms like “human”, “participants”, “infants” and includes numeric data (e.g., dosage or sample size).
	•	Animal-related: Mentions animals (e.g., “rat”, “mouse”, “rabbit”) and includes dosage information (e.g., “5 mg”, “100 µg”).
```
⸻

# 📦 Output
```
After execution, the script generates the following Excel files in the Scrapper_results folder:

File	Description
Scraped_Articles_Trimester.xlsx	 -> Raw scraped data for all keywords
Curated_Compounds_Trimester.xlsx ->	Filtered articles with experimental context
```

⸻