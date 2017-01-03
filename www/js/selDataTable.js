window.doNotChangeMyVar=false;

function getExploreTable () {
  return $('#exploreTable table').dataTable({'bRetrieve':true});
};

// Set the classes that TableTools uses to something suitable for Bootstrap
$.extend( true, $.fn.DataTable.TableTools.classes, {
  "container": "btn-group",
	"buttons": {
		"normal": "btn",
		"disabled": "btn disabled"
	},
	"collection": {
		"container": "DTTT_dropdown dropdown-menu",
		"buttons": {
			"normal": "",
			"disabled": "disabled"
		}
	}
} );

// Have the collection use a bootstrap compatible dropdown
$.extend( true, $.fn.DataTable.TableTools.DEFAULTS.oTags, {
	"collection": {
		"container": "ul",
		"button": "li",
		"liner": "a"
	}
} );

$(document).on('click', 'div.selectable table tbody tr', function(){
  var el = $(this);
  if((el.has('td.details').length > 0) || (window.doNotChangeMyVar)) {return}
  el.siblings().removeClass("rowsSelected");
  el.addClass("rowsSelected", this.clicked);
  el.trigger("change");
  var out = [];
  el.find("td").each(function(cell){
    out[cell] = $(this).text();
  });
  var ts = new Date().getTime();
//  var row = el.index()+1;
  var oTable = getExploreTable();
  Shiny.onInputChange("exploreTableVariation", out[2]+'\t'+out[3]+'\t'+out[5]+'\t'+oTable.fnPagingInfo().iStart+'\t'+ts);
});

$(document).on('mouseover', 'div.selectable table tbody tr',
    function(){el = $(this)
               if (el.has('td.details').length > 0)  {return}
               el.addClass("rowsHover");
              }
);

$(document).on('mouseout', 'div.selectable table tbody tr',
    function(){el = $(this)
               if (el.has('td.details').length > 0) {return}
               el.removeClass("rowsHover")
              }
);

function fnFormatDetails ( nTr )
{
    var oTable = getExploreTable();
    var aData = oTable.fnGetData( nTr );
    var aNotes = new Array;
    if (aData[8] != '') {aNotes.push(aData[8])};
    if (aData[7] != '') {aNotes = aNotes.concat(aData[7].split('\n'))};
    var sOut = '<div class="row-fluid">'
              +'<div class="span1"><em><b>Name</b></em></div><div class="span2"><em><b>Curated class</b></em></div>'
              +'<div class="span2"><em><b>Date</b></em></div><div class="span7"><em><b>Notes</b></em></div></div>'
              +aNotes.map(function(x){var aFields=x.split("\t");
                                      return('<div class="row-fluid">'
                                            +aFields.map(function(y,i){var sp=1;
                                                                       if(i==1){sp=2};
                                                                       if(i==2){sp=2};
                                                                       if(i==3){sp=7};
                                                                       y1=y.replace(/PMID:(\s+)?(\d+)/g,'<A HREF="http://www.ncbi.nlm.nih.gov/pubmed?term=$2" target="_blank">PMID:$2</A>');
                                                                       var aLines=y1.split("<br>");
                                                                       return('<div class="span'+sp+'">'
                                                                             +aLines.map(function(z){return('<div>'+z+'</div>')}).join('')
                                                                             +'</div>')
                                                                      }
                                                        ).join('')
                                            +'</div>');
                                     }
                          ).join('');
    return sOut;
};

$(document).on('mouseover', '#exploreTable tbody td img', function () {
  window.doNotChangeMyVar = true;
});

$(document).on('mouseout', '#exploreTable tbody td img', function () {
  window.doNotChangeMyVar = false;
});
  
$(document).on('click', '#exploreTable tbody td img', function () {
  var nTr = $(this).parents('tr')[0];
  var oTable = getExploreTable()
  if ( oTable.fnIsOpen(nTr) )
    {
      this.src = "images/details_open.png";
      oTable.fnClose( nTr );
    }
  else
    {
      this.src = "images/details_close.png";
      oTable.fnOpen( nTr, fnFormatDetails(nTr), 'details' );
    }
});

// Add a show all and hide all notes button to table tools
$.extend( true, $.fn.dataTable.TableTools.DEFAULTS, {
                  "sSwfPath": "TableTools-2.2.0/swf/copy_csv_xls_pdf.swf",
                  "aButtons": [{"sExtends":"text",
                                "sButtonText":"Show all notes",
                                 "fnClick": function(button,config){
                                              var oTable=getExploreTable();
                                              $('#exploreTable tbody tr td img[src$="open.png"]').each(function(){
                                                var nTr=$(this).parents('tr')[0];
                                                this.src='images/details_close.png';
                                                oTable.fnOpen(nTr,fnFormatDetails(nTr),'details');
                                                });
                                            }
                                },
                                {"sExtends":"text",
                                 "sButtonText":"Hide all notes",
                                 "fnClick": function(button,config){
                                              var oTable=getExploreTable();
                                              $('#exploreTable tbody tr td img[src$="close.png"]').each(function() {
                                                var nTr=$(this).parents('tr')[0];
                                                this.src='images/details_open.png';
                                                oTable.fnClose(nTr);
                                                });
                                            }
                                }]
                  }
);

//get rid of the CSV,PDF,Print
$.fn.DataTable.TableTools.DEFAULTS.aButtons.length=2;

/*
$.fn.dataTableExt.oApi.markRow = function ( oSettings, iRow, sClass )
{
  var el=$('tbody tr:nth-child('+iRow+')',oSettings.nTable);
  el.addClass(sClass);
}
*/
