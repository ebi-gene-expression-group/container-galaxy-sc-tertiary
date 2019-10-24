# Generation of tables

For generating the S2 table, you will need Perl XML:XPath, for instance you could use conda:

```bash
conda create -n perl-xpath perl-xml-xpath
```

Then activate the environment

```bash
conda activate perl-xpath
```

and create the table:

```bash
bash create_table_S2.sh
```
