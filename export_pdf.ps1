# PowerShell script to convert DOCX to PDF via Word COM Automation
$docxPath = "c:\Users\a3188798\OneDrive - Adelaide University\Desktop\NS BA AG\Bile_Acid_Analysis_Report.docx"
$pdfPath = "c:\Users\a3188798\OneDrive - Adelaide University\Desktop\NS BA AG\Bile_Acid_Analysis_Report.pdf"

try {
    # Remove existing PDF if present
    if (Test-Path $pdfPath) {
        Remove-Item $pdfPath -Force
    }

    Write-Output "Starting Microsoft Word..."
    $word = New-Object -ComObject Word.Application
    $word.Visible = $false
    $word.DisplayAlerts = 0

    Write-Output "Opening document..."
    # Open NOT as ReadOnly so ExportAsFixedFormat works
    $doc = $word.Documents.Open($docxPath, $false, $false)

    Write-Output "Exporting to PDF..."
    # wdExportFormatPDF = 17
    $doc.ExportAsFixedFormat($pdfPath, 17)

    Write-Output "Closing..."
    $doc.Close(0)
    $word.Quit()

    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($word) | Out-Null

    Write-Output "PDF export completed successfully!"
} catch {
    Write-Error "Error: $_"
    try { if ($doc) { $doc.Close(0) } } catch {}
    try { if ($word) { $word.Quit() } } catch {}
}
