#!/usr/bin/env python3
"""
selen_scraper.py - Selenium Edge Web Scraper for UQFF Data Collection
======================================================================

Scrapes astronomical data from websites without APIs using Edge WebDriver.

Usage:
    python selen_scraper.py --url "https://example.com" --output data.csv
    python selen_scraper.py --config scraper_config.json

Author: Daniel T. Murphy
Framework: UQFF Star-Magic
"""

import os
import sys
import time
import json
import csv
import argparse
from datetime import datetime
from typing import List, Dict, Any, Optional
from pathlib import Path

from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException


class SeleniumScraper:
    """
    Selenium Edge WebDriver scraper for astronomical data collection.
    """
    
    def __init__(self, headless: bool = True, timeout: int = 30):
        """
        Initialize the Edge WebDriver.
        
        Args:
            headless: Run browser without GUI
            timeout: Default wait timeout in seconds
        """
        self.timeout = timeout
        self.driver = None
        
        # Edge driver path
        edge_driver_path = os.path.join(os.getcwd(), 'edge_driver', 'msedgedriver.exe')
        
        # Edge options
        edge_options = Options()
        if headless:
            edge_options.add_argument('--headless')
        edge_options.add_argument('--disable-gpu')
        edge_options.add_argument('--no-sandbox')
        edge_options.add_argument('--disable-dev-shm-usage')
        edge_options.add_argument('--window-size=1920,1080')
        edge_options.add_argument('--disable-blink-features=AutomationControlled')
        edge_options.add_experimental_option("excludeSwitches", ["enable-automation"])
        edge_options.add_experimental_option('useAutomationExtension', False)
        
        # Initialize Edge service and driver
        service = Service(executable_path=edge_driver_path)
        self.driver = webdriver.Edge(service=service, options=edge_options)
        self.driver.execute_script("Object.defineProperty(navigator, 'webdriver', {get: () => undefined})")
        
        print(f"✓ Edge WebDriver initialized (headless={headless})")
    
    def get_page(self, url: str, wait_for: Optional[str] = None) -> bool:
        """
        Navigate to URL and optionally wait for element.
        
        Args:
            url: Target URL
            wait_for: CSS selector to wait for (optional)
            
        Returns:
            True if successful, False otherwise
        """
        try:
            print(f"→ Loading: {url}")
            self.driver.get(url)
            
            if wait_for:
                WebDriverWait(self.driver, self.timeout).until(
                    EC.presence_of_element_located((By.CSS_SELECTOR, wait_for))
                )
                print(f"✓ Element found: {wait_for}")
            else:
                time.sleep(2)  # Brief wait for page load
            
            return True
        except TimeoutException:
            print(f"✗ Timeout waiting for element: {wait_for}")
            return False
        except Exception as e:
            print(f"✗ Error loading page: {e}")
            return False
    
    def extract_table(self, table_selector: str = "table") -> List[List[str]]:
        """
        Extract data from HTML table.
        
        Args:
            table_selector: CSS selector for table element
            
        Returns:
            List of rows (each row is a list of cell values)
        """
        try:
            table = self.driver.find_element(By.CSS_SELECTOR, table_selector)
            rows = table.find_elements(By.TAG_NAME, "tr")
            
            data = []
            for row in rows:
                cells = row.find_elements(By.TAG_NAME, "td")
                if not cells:  # Try th for headers
                    cells = row.find_elements(By.TAG_NAME, "th")
                
                row_data = [cell.text.strip() for cell in cells]
                if row_data:  # Skip empty rows
                    data.append(row_data)
            
            print(f"✓ Extracted {len(data)} rows from table")
            return data
        except NoSuchElementException:
            print(f"✗ Table not found: {table_selector}")
            return []
        except Exception as e:
            print(f"✗ Error extracting table: {e}")
            return []
    
    def extract_elements(self, selector: str, attribute: Optional[str] = None) -> List[str]:
        """
        Extract text or attribute from elements matching selector.
        
        Args:
            selector: CSS selector for elements
            attribute: Attribute name to extract (None for text content)
            
        Returns:
            List of extracted values
        """
        try:
            elements = self.driver.find_elements(By.CSS_SELECTOR, selector)
            
            if attribute:
                data = [elem.get_attribute(attribute) for elem in elements]
            else:
                data = [elem.text.strip() for elem in elements]
            
            data = [d for d in data if d]  # Remove empty strings
            print(f"✓ Extracted {len(data)} values from '{selector}'")
            return data
        except Exception as e:
            print(f"✗ Error extracting elements: {e}")
            return []
    
    def click_element(self, selector: str) -> bool:
        """
        Click element matching selector.
        
        Args:
            selector: CSS selector for element to click
            
        Returns:
            True if successful, False otherwise
        """
        try:
            element = WebDriverWait(self.driver, self.timeout).until(
                EC.element_to_be_clickable((By.CSS_SELECTOR, selector))
            )
            element.click()
            print(f"✓ Clicked: {selector}")
            time.sleep(1)  # Brief wait after click
            return True
        except Exception as e:
            print(f"✗ Error clicking element: {e}")
            return False
    
    def input_text(self, selector: str, text: str, submit: bool = False) -> bool:
        """
        Input text into form field.
        
        Args:
            selector: CSS selector for input field
            text: Text to input
            submit: Submit form after input
            
        Returns:
            True if successful, False otherwise
        """
        try:
            element = WebDriverWait(self.driver, self.timeout).until(
                EC.presence_of_element_located((By.CSS_SELECTOR, selector))
            )
            element.clear()
            element.send_keys(text)
            print(f"✓ Input text into '{selector}': {text}")
            
            if submit:
                element.submit()
                print(f"✓ Form submitted")
                time.sleep(2)  # Wait for results
            
            return True
        except Exception as e:
            print(f"✗ Error inputting text: {e}")
            return False
    
    def get_page_source(self) -> str:
        """Get current page HTML source."""
        return self.driver.page_source
    
    def screenshot(self, filename: str) -> bool:
        """
        Save screenshot of current page.
        
        Args:
            filename: Output filename
            
        Returns:
            True if successful, False otherwise
        """
        try:
            self.driver.save_screenshot(filename)
            print(f"✓ Screenshot saved: {filename}")
            return True
        except Exception as e:
            print(f"✗ Error saving screenshot: {e}")
            return False
    
    def close(self):
        """Close the browser and quit the driver."""
        if self.driver:
            self.driver.quit()
            print("✓ Browser closed")


def save_to_csv(data: List[List[str]], filename: str):
    """Save scraped data to CSV file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"{filename}_{timestamp}.csv"
    
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerows(data)
    
    print(f"✓ Data saved to: {output_file}")
    return output_file


def save_to_json(data: Any, filename: str):
    """Save scraped data to JSON file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"{filename}_{timestamp}.json"
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print(f"✓ Data saved to: {output_file}")
    return output_file


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE SCRAPING TASKS
# ═══════════════════════════════════════════════════════════════════════════════

def scrape_arxiv_papers(query: str, max_results: int = 50) -> List[Dict[str, str]]:
    """
    Scrape arXiv for papers matching query.
    
    Args:
        query: Search query
        max_results: Maximum number of results
        
    Returns:
        List of paper metadata dictionaries
    """
    scraper = SeleniumScraper(headless=False)
    papers = []
    
    try:
        # Navigate to arXiv
        url = f"https://arxiv.org/search/?query={query}&searchtype=all&abstracts=show&size={max_results}"
        scraper.get_page(url, wait_for=".arxiv-result")
        
        # Extract paper data
        titles = scraper.extract_elements(".title.is-5", attribute=None)
        authors = scraper.extract_elements(".authors", attribute=None)
        abstracts = scraper.extract_elements(".abstract-full", attribute=None)
        
        for i in range(min(len(titles), max_results)):
            papers.append({
                'title': titles[i] if i < len(titles) else '',
                'authors': authors[i] if i < len(authors) else '',
                'abstract': abstracts[i] if i < len(abstracts) else ''
            })
        
        print(f"✓ Scraped {len(papers)} papers from arXiv")
        
    finally:
        scraper.close()
    
    return papers


def scrape_astronomical_catalog(url: str, table_selector: str = "table") -> List[List[str]]:
    """
    Scrape astronomical data from catalog website.
    
    Args:
        url: Catalog URL
        table_selector: CSS selector for data table
        
    Returns:
        Scraped table data
    """
    scraper = SeleniumScraper(headless=False)
    data = []
    
    try:
        scraper.get_page(url, wait_for=table_selector)
        data = scraper.extract_table(table_selector)
        print(f"✓ Scraped {len(data)} rows from catalog")
        
    finally:
        scraper.close()
    
    return data


# ═══════════════════════════════════════════════════════════════════════════════
# CLI INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='Selenium Edge Web Scraper for UQFF')
    parser.add_argument('--url', type=str, help='Target URL to scrape')
    parser.add_argument('--table', type=str, default='table', help='CSS selector for table')
    parser.add_argument('--output', type=str, default='scraped_data', help='Output filename')
    parser.add_argument('--arxiv', type=str, help='Search arXiv for papers')
    parser.add_argument('--headless', action='store_true', help='Run in headless mode')
    parser.add_argument('--max-results', type=int, default=50, help='Maximum results')
    
    args = parser.parse_args()
    
    if args.arxiv:
        # Scrape arXiv
        papers = scrape_arxiv_papers(args.arxiv, args.max_results)
        save_to_json(papers, args.output)
    
    elif args.url:
        # Scrape table from URL
        data = scrape_astronomical_catalog(args.url, args.table)
        if data:
            save_to_csv(data, args.output)
    
    else:
        # Interactive mode
        print("\n" + "="*80)
        print("SELENIUM EDGE SCRAPER - UQFF STAR-MAGIC")
        print("="*80)
        print("\nExamples:")
        print("  1. Scrape arXiv papers:")
        print("     python selen_scraper.py --arxiv 'black hole UQFF' --max-results 20")
        print("\n  2. Scrape data table:")
        print("     python selen_scraper.py --url 'https://example.com/data' --table 'table.data'")
        print("\n  3. Custom scraping (modify script for your needs)")
        print("="*80 + "\n")


if __name__ == "__main__":
    main()
