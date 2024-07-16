#Long Read Alignment Plotting

To generate plots, run a command with the following form

```bash
python3 ./import/prepare-alignments.py -L /path/to/my.bam chr1:91665905-91667512 > data.json
```

Then in the `axaplotl` directory, run

```bash
npm install
npm run dev
```

Which will create an app running on port 3000, which you can point your browser to (i.e. <http://localhost:3000/>). In the app, use the file picker to select the json file you created with the alignments, and watch the plots appear!
