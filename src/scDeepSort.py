# run scDeepSort on human blood data

import deepsort
z = deepsort.DeepSortPredictor("human", "Blood")
print(z.predict("../data/scDeepSort/human_Blood3223_data.csv"))
