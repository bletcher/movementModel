# movementModel

Movement models for PIT tagged fish in streams

Details:
1. Model runs take a long time with all observations including those after the last observation (full dataset).
2. Limit movement model parameter estimation to knownZ == 1 and to fish with more than 1 observation (reduced dataset).
3. Simulate movement for the full dataset (knownZ == [0,1]) based on parameter estimates from the reduced dataset.
4. Use the simulated locations as input to the growth model (https://github.com/bletcher/growthModel).


Original Steps: get and prepare raw data, viz raw data, iterate to get antenna data correct (in WB), develop model for river transitions including captures, antenna data, wanding data, and acoustic tag data. Start with SB then do WB. After river-based model, do a section-within-river model. This will add habitat and density covariates. Finally feed the movement model output into the growth model, then feed growth model results into survival model.
