selDataTableOutput <- function (outputId) 
{
  tagList(
    singleton(
      tags$head(
        tags$link(href="shared/datatables/css/DT_bootstrap.css", rel="stylesheet", type="text/css"),
        tags$link(href="TableTools-2.2.0/css/dataTables.tableTools.min.css", rel="stylesheet", type="text/css"),
        tags$style(type="text/css",".rowsHover td{background-color: rgba(245, 222, 179, 0.4) !important}
                                    .visited td{background-color: rgba(232, 101, 55, 0.6) !important}
                                    .visitedLast td{background-color: rgba(232, 101, 55, 0.8) !important}
                                    .selectable div table tbody tr{cursor: hand; cursor: pointer;}
                                    .details {background-color: #D1CFD0 !important;border: 1px solid #A19B9E;}
                                    div.dataTables_length label {float: right !important; text-align: right !important;}
                                    div.dataTables_info {float: right ! important;}
                                    div.dataTables_filter label {float: left !important;}
                                    div.dataTables_paginate {float: left !important; margin: 0;}"
        ),
        HTML('<style type="text/css">                                
              .pagination ul>li>a, .pagination ul>li>span, a.btn>span, .dataTables_info {
                 -webkit-user-select: none;}
              </style>'),
        tags$script(src = "shared/datatables/js/jquery.dataTables.min.js"),
        tags$script(src = "TableTools-2.2.0/js/dataTables.tableTools.min.js"),
        tags$script(src = "shared/datatables/js/DT_bootstrap.js"),
        tags$script(src = "js/selDataTable.js")
      )
    ), 
    div(id = outputId, class = "shiny-datatable-output selectable")
  )
}
