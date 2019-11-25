import grakn
from tabulate import tabulate

client = grakn.Grakn(uri = "localhost:48555")

with client.session(keyspace = "biograkn_dn") as session:
    
    with session.transaction(grakn.TxType.READ) as tx:
        
        uniprot_id_type = tx.get_schema_concept("uniprot-id")
       
        print("\n\n+ Protein(s) encoded by the gene with entrez-id of 100137049:\n", flush=True)

        it = tx.query('match $gpe (encoding-gene: $ge, encoded-protein: $pr) isa gene-protein-encoding; $ge isa gene has entrez-id "100137049"; limit 10; get;', infer=False)
        
        prots = list()
        
        for res in it:
            prots.append([next(res.map()['pr'].attributes(uniprot_id_type)).value()])
        
        print(tabulate(prots, headers=["Protein"]))

        
    with session.transaction(grakn.TxType.WRITE) as tx:
        
        print("\n\n+ Diseases affecting the appendix tissue:\n", flush=True)
        
        it = tx.query('match $ti isa tissue has tissue-name "appendix"; $dta (associated-disease: $di, associated-tissue: $ti) isa disease-tissue-association; $di isa disease; limit 10; get;')
        
        dis = list()
        for res in it:
            dis.append([next(res.map()['di'].attributes()).value()])
        print(tabulate(dis, headers=["Disease"]))


        print("\n\n+ Proteins associated with Asthma disease:\n", flush=True)
        
        it = tx.query('match $di isa disease has disease-name "Asthma"; $dda (associated-protein: $pr, associated-disease: $di) isa protein-disease-association; limit 10; get;', infer=False)
        
        prots = list()
        
        for res in it:
            prots.append([next(res.map()['pr'].attributes(uniprot_id_type)).value()])
        
        print(tabulate(prots, headers=["Protein"]))


    with session.transaction(grakn.TxType.READ) as tx:

        print("\n\n+ Diseases associated with protein interactions taking place in the liver:\n", flush=True)
        
        uniprot_id_type = tx.get_schema_concept("uniprot-id")

        it = tx.query('match $ti isa tissue, has tissue-name "liver"; $di isa disease; $pl (tissue-context: $ti, biomolecular-process: $ppi) isa process-localisation; $ppi (interacting-protein: $pr, interacting-protein: $pr2) isa protein-protein-interaction; $pr isa protein; $pr2 isa protein; $pr != $pr2; $pda (associated-protein: $pr, associated-disease: $di) isa protein-disease-association; limit 30; get;', infer=False)
        
        table = list()
        
        for res in it:
            table.append([next(res.map()['pr'].attributes(uniprot_id_type)).value() + " <-> " + next(res.map()['pr2'].attributes(uniprot_id_type)).value(), next(res.map()['di'].attributes()).value()])

        print(tabulate(table, headers=["Protein-Protein interaction","Disease"], colalign=("center",)))


    with session.transaction(grakn.TxType.WRITE) as tx:

        print("\n\n+ Drugs and diseases associated with the same differentially expressed gene from comparisons made in geo-series with id of GSE27876:\n", flush=True)
        
        gene_symbol_type = tx.get_schema_concept("gene-symbol")
        drug_name_type = tx.get_schema_concept("drug-name")
        
        it = tx.query('match $geo-se isa geo-series has GEOStudy-id "GSE27876"; $comp (compared-groups: $geo-comp, containing-study: $geo-se) isa comparison; $def (conducted-analysis: $geo-comp, differentially-expressed-gene: $ge) isa differentially-expressed-finding; $dgi (target-gene: $ge, interacted-drug: $dr) isa drug-gene-interaction; $gda (associated-gene: $ge, associated-disease: $di) isa gene-disease-association; limit 10; get;')
        
        table = list()
        
        for res in it:
            table.append([next(res.map()['ge'].attributes(gene_symbol_type)).value(), next(res.map()['dr'].attributes(drug_name_type)).value(), next(res.map()['di'].attributes()).value()])

        print(tabulate(table, headers=["Gene","Drug","Disease"], colalign=("center",)))
