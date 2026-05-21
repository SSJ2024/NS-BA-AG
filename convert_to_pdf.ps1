# PowerShell script to convert DOCX to PDF via Word COM Automation
$docxPath = "c:\Users\a3188798\OneDrive - Adelaide University\Desktop\NS BA AG\Bile_Acid_Analysis_Report.docx"
$pdfPath = "c:\Users\a3188798\OneDrive - Adelaide University\Desktop\NS BA AG\Bile_Acid_Analysis_Report.pdf"

try {
    Write-Output "Starting Microsoft Word application in background..."
    $word = New-Object -ComObject Word.Application
    $word.Visible = $false
    $word.DisplayAlerts = 0 # wdAlertsNone = 0

    Write-Output "Opening document: $docxPath"
    # Open(FileName, ConfirmConversions, ReadOnly)
    $doc = $word.Documents.Open($docxPath, $false, $true)

    Write-Output "Exporting to PDF: $pdfPath"
    # SaveAs(FileName, FileFormat)
    # 17 represents wdFormatPDF
    $doc.SaveAs($pdfPath, 17)

    Write-Output "Closing document..."
    $doc.Close(0) # wdDoNotSaveChanges = 0

    Write-Output "Quitting Microsoft Word..."
    $word.Quit()

    Write-Output "PDF conversion completed successfully!"
} catch {
    Write-Error "An error occurred during PDF conversion: $_"
}
