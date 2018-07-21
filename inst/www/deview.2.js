Ext.define('MyStore', {
    extend: 'Ext.data.Store',
    autoLoad: true,
    fields: ['gene','membrane', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj','Z','Za','pvalue','padj','significant'],
    pageSize: 50,
    remoteFilter: true,
    remotePaging: true,
    remoteSort: true
})

// clear filter filed trigger button
Ext.define('Ext.ux.CustomTrigger', {
    extend: 'Ext.form.field.Trigger',
    alias: 'widget.customtrigger',
    initComponent: function () {
        var me = this;
        me.triggerCls = 'x-form-clear-trigger';
        me.callParent(arguments);
    },
    // override onTriggerClick
    onTriggerClick: function() {
        if(this.getValue()!='') {
            this.setRawValue('');
            this.fireEvent('change',this,'');
        }
    }
});

var getUrlParameter = function getUrlParameter(sParam) {
    var sPageURL = decodeURIComponent(window.location.search.substring(1)),
        sURLVariables = sPageURL.split('&'),
        sParameterName,
        i;

    for (i = 0; i < sURLVariables.length; i++) {
        sParameterName = sURLVariables[i].split('=');

        if (sParameterName[0] === sParam) {
            return sParameterName[1] === undefined ? true : sParameterName[1];
        }
    }
};

var store = Ext.create('MyStore')

var ptbar = Ext.create('Ext.PagingToolbar', {
        store: store,
        displayInfo: true,
	
        //displayMsg: 'Displaying genes {0} - {1} of {2}',
        emptyMsg: "No DE genes to display",
        items:[
            // {
            //     flex:5,
            //     width: 100,
            //     minWidth: 50,
            //     xtype: 'customtrigger',
            //     emptyText: 'filter by gene name...',
            //     listeners: {
            //         change: {buffer: 200, fn: function(field, value) {
            //             if (value.length>0) {
            //                 store.clearFilter(true);
            //                 store.filter({property: 'gene', value: value});

	    // 			store.load();
            //             } else {
            //                 store.clearFilter(false);
	    // 		    store.load();
            //             }
            //         }}
            //     }
            //},
	    '-',
	    {
		fieldLabel: 'Page size',
		flex:1,
		xtype: 'numberfield',
	        tooltip: 'Number of genes to show per page',
		label: 'N genes',
		minWidth: 80,
		width: 80,
		value: 50,
		minValue: 1,
		maxValue: 10000,
		disabled: false,
		listeners : {
		    change : {buffer: 800, fn:function(f,v) {
			store.pageSize=parseInt(v);
			store.load();
		    }}
		}
	    }
	],
        listeners: {
            afterrender: function() {
                this.down('#refresh').hide();
            }
        }
})


function OpenWindowWithPost(url, windowoption, name, params)
{
    var form = document.createElement("form");
    form.setAttribute("method", "post");
    form.setAttribute("action", url);
    form.setAttribute("target", name);
    form.setAttribute("enctype", "multipart/form-data");
    
    for (var i in params) {
	var input = document.createElement('input');
	input.type = 'hidden';
	input.name = params[i].Key;
	input.value = params[i].Value;
	form.appendChild(input);
    }
    
    document.body.appendChild(form);
    
    //note I am using a post.htm page since I did not want to make double request to the page 
    //it might have some Page_Load call which might screw things up.
    window.open("post.htm", name, windowoption);
    
    form.submit();
    
    document.body.removeChild(form);
}

Ext.require('Ext.ux.grid.FiltersFeature');

Ext.define('MyGrid', {
    extend: 'Ext.grid.Panel',
    //title: 'Differentially Expressed Genes',
    header: false,
    width: 800,
    height: 500,
    features: [{
        ftype: 'filters',
	local: false
    }],
    
    // show boxplots upon selection
    selType: 'rowmodel',
    listeners: { 
	selectionchange: function(selModel,selected) {
	    if(selected.length>0) {
		var rid=selected[0].raw['rowid'];
		var data=[];
		Object.keys(dedata.ilev).forEach(function(key) {
		    var res = {
			text: dedata.ilev[key].snames,
			y: dedata.ilev[key].val[rid-1],
			name: key,
			type: 'box',
			jitter: 0.3,
			pointpos: -1.8,
			boxpoints: 'all'
		    }
		    data.push(res);
		})
		Plotly.newPlot('plots', data,plotLayout);
	    } else {
		console.log("boo")
	    }
	    
	}
    },

    columns: [{
        text: 'gene',
        dataIndex: 'gene',
	flex: 1,
	renderer: function(value,m,r) {
	    return Ext.String.format('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene={0}" target="_blank">{1}</a>',value,value)
	},
	filter: { type: 'string' },
    // }, {
    //     text: 'baseMean',
    //     dataIndex: 'baseMean',
    // 	xtype: 'numbercolumn',
	
    // 	format: '0.00'
    }, {
	text: 'membrane',
	dataIndex: 'membrane',
	xtype: 'booleancolumn',
	filter: { type: 'boolean' }
	//filter: { type: 'string' }
    }, {
        text: 'log2 Fold Change',
        dataIndex: 'log2FoldChange',
	xtype: 'numbercolumn',
	filter: { type: 'numeric' },
	format: '0.00'
    }, {
        text: 'Z',
        dataIndex: 'Z',
	xtype: 'numbercolumn',
	format: '0.00'
    }, {
        text: 'Zadj',
        dataIndex: 'Za',
	xtype: 'numbercolumn',
	filter: { type: 'numeric' },
	format: '0.00'
    }, /*{
	text: 'P-value',
	dataIndex: 'pvalue',
	xtype: 'numbercolumn',
	filter: { type: 'numeric' }
    }, {
	text: 'Adj P-value',
	dataIndex: 'padj',
	xtype: 'numbercolumn',
	filter: { type: 'numeric' }
    },*/ {
	text: 'Significant',
	dataIndex: 'significant',
	xtype: 'booleancolumn',
	filter: { type: 'boolean' }
    }




],
    tbar: ptbar,
    bbar: [
	{
	    xtype: 'button',
	    text: 'Run Enrichr',
	    id: 'te',
	    tooltip: 'Test enrichment of the genes on the current page using Enrichr',
	    // var records = store.getRange();`
	    handler: function(b,e) {
		var records=store.getRange();
		var genes=new Array(records.length);
		for(var i=0; i<records.length; i++ ) {
		    genes[i]=records[i].get('gene');
		}
		
		//window.open("http://amp.pharm.mssm.edu/Enrichr/enrich?list="+genes.join("\n"), "enrichment", "");
		var param = [{ Key: 'list', Value: genes.join("\n")}, { Key : 'species', Value: 'human_hg19'}];
		OpenWindowWithPost("http://amp.pharm.mssm.edu/Enrichr/enrich","","enrich",param);

		//console.log(genes.join("\n"));
	    }
	},
	{
	    xtype: 'button',
	    text: 'Run GOrilla',
	    id: 'teg',
	    tooltip: 'Test enrichment of the genes on the current page using GOrilla',
	    // var records = store.getRange();`
	    handler: function(b,e) {
		var records=store.getRange();
		var genes=new Array(records.length);
		for(var i=0; i<records.length; i++ ) {
		    genes[i]=records[i].get('gene');
		}
		
		//window.open("http://amp.pharm.mssm.edu/Enrichr/enrich?list="+genes.join("\n"), "enrichment", "");
		var param = [{ Key: 'target_set', Value: genes.join("\n")}, {Key: 'application', Value: 'gorilla'}, { Key : 'species', Value: 'HOMO_SAPIENS'}, {Key: 'run_mode', Value: 'hg'}, {Key: 'background_set', Value: dedata.genes.join("\n")}, {Key: 'db', Value: 'proc'}, {Key: 'pvalue_thresh', Value: '0.001'}];
		OpenWindowWithPost("http://cbl-gorilla.cs.technion.ac.il/servlet/GOrilla","","enrich",param);

		//console.log(genes.join("\n"));
		//console.log("background set:")
		//console.log(dedata.genes.join("\n"));
	    }
	}
    ]

})
var dedata;

var plotLayout = {
    yaxis: {
	autorange: true,
	showgrid: true,
	zeroline: true,
	title: 'log(expression)',
    },
    margin: {
	l: 50,
	r: 50,
	b: 50,
	t: 50
    },
    showlegend: false
};

Ext.application({
    name: 'Fiddle',
    
    launch: function() {
	
	// load data
	Ext.Ajax.request({
	    url: getUrlParameter('d'),
	    success: function(response, opts) {
		dedata = Ext.decode(response.responseText);
		console.log("boo")
		store.setProxy({type: "memory",enablePaging: true, data: dedata.res})
		store.load();
	    },
	    failure: function(response, opts) {
		console.log('server-side failure with status code ' + response.status);
	    }
	});
        
	
        var myGridPanel = Ext.create('MyGrid', {
            store: store
        })


	var viewport = Ext.create('Ext.Viewport', {
            layout: {
		type: 'border',
		padding: 5
            },
            defaults: {
		split: true
            },
            items: [{ region: 'west',
		      width: 800,
		      layout: 'fit',
		      items: [myGridPanel]
		    },
		    { region: 'center',
		      xtype: 'panel',
		      html:['<div id="plots"/>' ]
		    }]
	}).show();

    }
});
