import sys
import pdfkit

def generate_pdf(html_file, output_pdf):
    try:
        # Configuration : assure-toi que wkhtmltopdf est bien installé
        config = pdfkit.configuration(wkhtmltopdf='/usr/bin/wkhtmltopdf')  # chemin Linux par défaut
        pdfkit.from_file(html_file, output_pdf, configuration=config)
        print(f"[✔] PDF report generated: {output_pdf}")
    except Exception as e:
        print(f"[✘] PDF generation failed: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_pdf_report.py <HTML_FILE> <OUTPUT_PDF>")
        sys.exit(1)

    html_input = sys.argv[1]
    pdf_output = sys.argv[2]

    generate_pdf(html_input, pdf_output)
