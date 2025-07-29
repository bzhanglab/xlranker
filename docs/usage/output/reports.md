# XLRanker Reports

XLRanker generates a variety of reports at different levels of confidence.

1. Conservative
2. Minimal
3. Expanded
4. All-inclusive

## Conservative Report

The conservative report is the most stringent and includes only the highest confidence interactions. These interactions are protein pairs that are part of a single-membered protein group that was parsimoniously selected. _This level may result in not all peptide sequences being represented._

## Minimal Report

The minimal report includes all the interactions in the minimal report as well as the primarily selected interactions selected from the machine learning model, as well as best intra pair from parsimoniously selected protein groups. This will level will have **all** peptide sequences represented.

## Expanded Report

This report includes all interactions selected during the parsimonious and machine learning selection steps.