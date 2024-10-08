function showLoading() {
    document.getElementById('loader').style.display = 'block';
}

function checkGEXAndSubmit() {
    var gexRadio = document.getElementById('GEX');
    if (gexRadio.checked) {
        alert('Uploading only GEX data is not supported yet.');
        return false; // Prevent form submission
    }
    return true; // Allow form submission
}