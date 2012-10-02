
samples = set(['tt','wj','mj','st','dy'])

luminosity = (5008, 0.04) # (1/pb, %)

xs = {'tt' : ( 149.600, 0.40), # (pb, %)
      'wj' : (1911.800, 0.40),
      'mj' : (   1.000, None),
      'st' : (  71.968, 0.01),
      'dy' : (2475.000, 0.01)}

components = dict( [(item,('',1.0)) for item in ['wj','mj','st','dy']] +
                   [('tt',[('gg',0.80),
                           ('qg',0.10),
                           ('qq',0.07),
                           ('ag',0.03)])] )

efficiency = {'mu' : { 'wj' : [('',1.0)],
                       'mj' : [('',1.0)],
                       'st' : [('',1.0)],
                       'dy' : [('',1.0)],
                       'tt' : [('gg',1.0),
                               ('qg',1.0),
                               ('qq',1.0),
                               ('ag',1.0)] 
                       },
              'el' : { 'wj' : [('',1.0)],
                       'mj' : [('',1.0)],
                       'st' : [('',1.0)],
                       'dy' : [('',1.0)],
                       'tt' : [('gg',1.0),
                               ('qg',1.0),
                               ('qq',1.0),
                               ('ag',1.0)] 
                       }
              }

def histogram(dist, chan, samp, comp = '') : 
    import random, ROOT as r
    dummy_hist = r.TH1D('dummy'+((4*'%s')%tuple([random.randrange(0,10) for i in range(4)])),'',100,0,1)
    sigma = random.random()
    for i in range(10000) : dummy_hist.Fill(random.gauss(0.5,sigma))
    return dummy_hist

