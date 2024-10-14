from bokeh.models import CustomJS

def side_panel_callback(source, x, project_id):
    base_url = "/analyze/gene_agg/"  # Base URL for the Flask route

    return CustomJS(args=dict(source=source, base_url=base_url, project_id=project_id), code="""
        var selected = source.selected.indices;
        if (selected.length > 0) {
            var index = selected[0];
            var gene = source.data['""" + x + """'][index];
            var count = source.data['counts'][index];
            var content = "<h3>Gene: " + gene + "</h3>" +
                        "<p>Count: " + count + "</p>";
                        
            var options = ['j_call_VDJ', 'v_call_VDJ', 'c_call_VDJ', 'v_call_VJ', 'j_call_VJ', 'c_call_VJ'];

            options.forEach(function(option) {
                if (option !== '""" + x + """') {
                    var data = source.data[option][index];
                    var entries = data.split('\\n');
                    content += "<p>" + option + ":</p>";
                    entries.forEach(function(entry) {
                        // Split the entry at the colon and take the first part
                        var gene_part = entry.split(':')[0].trim();
                        // Construct the URL directly in JavaScript
                        var url = base_url + encodeURIComponent(project_id) + '/' + encodeURIComponent(gene) + '/' + encodeURIComponent(gene_part);
                        content += "<pre><a href='" + url + "' target='_blank'>" + entry + "</a></pre><br>";
                    });
                }
            });

            openPanel(content);
        }
    """)