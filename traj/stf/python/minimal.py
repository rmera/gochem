def weaita(sele="test_stf",objname="vieja"):
    for i in range(2):
        b=cmd.get_model(sele)
        cmd.load_model(b,objname,state=i+1)




cmd.extend("weaita",weaita) 
