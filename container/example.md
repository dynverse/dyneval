Make sure docker is installed, and the folder ~/tmp exists (or change it in the code)

```
docker run -v ~/tmp/:/ti/ dynverse/dyntoy --model linear --num_cells 100 --output /ti/dataset.h5

docker run -v ~/tmp/:/ti/ dynverse/ti_slingshot --dataset /ti/dataset.h5 --output /ti/model.h5

docker run -v ~/tmp/:/ti/ dynverse/dyneval --dataset /ti/dataset.h5 --model /ti/model.h5 --output /ti/scores.json
```

The scores are in ~/tmp/scores.json
