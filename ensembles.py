

def ensemble_specs():
    base = {'lumiFactor':1.0, 'Nens':1000}

    ups = [{'label':'A','letter':'A'},
           {'label':'B','letter':'B'},
           {'label':'C','letter':'C'},
           {'label':'D','letter':'D'},
           {'label':'A4','letter':'A','lumiFactor':0.25, 'Nens':100},
           {'label':'A2','letter':'A','lumiFactor':0.50, 'Nens':100},
           {'label':'2A','letter':'A','lumiFactor':2.00, 'Nens':100},
           {'label':'4A','letter':'A','lumiFactor':4.00, 'Nens':100},
           {'label':'8A','letter':'A','lumiFactor':8.00, 'Nens':100},
    ]
    specs = [dict(base) for u in ups]
    for s,u in zip(specs,ups): s.update(u)
    return specs