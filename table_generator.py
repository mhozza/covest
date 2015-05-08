#! /usr/bin/env python


def format_table(header, lines, template, line_template, header_cell_template,
                 cell_template, round_floats=4):
    lbody = ''
    for c in header:
        lbody += header_cell_template.format(cbody=c)
    thead = line_template.format(lbody=lbody)
    tbody = ''
    for l in lines:
        lbody = ''
        for c in header:
            val = l.get(c, None)
            if round_floats and type(val) is float:
                val = round(val, round_floats)
            lbody += cell_template.format(cbody=val)
        tbody += line_template.format(lbody=lbody)
    return template.format(thead=thead, tbody=tbody)


def format_table_html(header, lines):
    template = '''<!DOCTYPE html>
    <html>
    <head>
    <!-- Latest compiled and minified CSS -->
    <link rel="stylesheet"
        href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css">
    </head>
    <table class="table table-condensed table-hover">
        {thead}
        {tbody}
    </table>
    </html>
    '''
    line_template = '\t<tr>\n{lbody}\t</tr>\n'
    cell_template = '\t\t<td>{cbody}</td>\n'
    header_cell_template = '\t\t<th>{cbody}</th>\n'
    return format_table(header, lines, template, line_template,
                        header_cell_template, cell_template)
