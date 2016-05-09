<!DOCTYPE html>
<html>
<head>
<!-- Latest compiled and minified CSS -->
<link rel="stylesheet"
  href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css">
</head>
<body>
<table class="table table-condensed table-hover">
    <thead>
        <tr>{{# header }}
            <th>{{ value }}</th>{{/ header }}
        </tr>
    </thead>
    <tbody>
        {{# body }}
        <tr>{{# line }}
            <td>{{ value }}</td>{{/ line }}
        </tr>
        {{/ body }}
    </tbody>
</table>
</body>
</html>
