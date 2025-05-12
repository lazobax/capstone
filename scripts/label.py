# # Interactive CONSEQ Labeling Script

# # List of unique CONSEQ values
# conseq_list = [
#     'Pathogenic', 'benign', 'Benign', 'Uncertain significance',
#        'benign, VUS', 'Likely benign',
#        'Conflicting classifications of pathogenicity', 'VUS',
#        'Benign/Likely benign', 'likely benign',
#        'benign, likely benign, VUS', '-', 'nan', 'likely benign, VUS',
#        'pathogenic, VUS', 'Likely pathogenic', 'benign, NA',
#        'pathogenic (dominant)', 'NA, VUS',
#        'pathogenic, pathogenic (dominant)', 'not provided',
#        'NA, pathogenic, pathogenic (dominant)', 'NA, pathogenic, VUS',
#        'Pathogenic/Likely pathogenic', 'pathogenic',
#        'pathogenic, pathogenic (dominant), VUS', 'likely benign, NA',
#        'likely pathogenic', 'likely pathogenic (dominant), NA',
#        'likely pathogenic, VUS', 'likely pathogenic, NA',
#        'likely pathogenic, NA, pathogenic, VUS',
#        'likely pathogenic (dominant), pathogenic, pathogenic (dominant)',
#        'benign, NA, VUS', 'likely pathogenic (dominant)',
#        'NA, pathogenic',
#        'likely pathogenic, NA, pathogenic, pathogenic (dominant)',
#        'pathogenic (dominant), VUS',
#        'likely pathogenic, likely pathogenic (dominant), NA, pathogenic, pathogenic (dominant), VUS',
#        'likely benign, likely pathogenic (dominant), pathogenic, VUS',
#        'benign, likely benign', 'benign, likely benign, NA, VUS',
#        'likely pathogenic, pathogenic, pathogenic (dominant)',
#        'likely pathogenic (dominant), pathogenic',
#        'likely benign, NA, VUS',
#        'likely pathogenic, pathogenic (dominant)',
#        'likely pathogenic, pathogenic',
#        'benign, NA, pathogenic, pathogenic (dominant), VUS',
#        'NA, pathogenic, pathogenic (dominant), VUS',
#        'benign, likely pathogenic', 'likely pathogenic, pathogenic, VUS',
#        'likely pathogenic (dominant), NA, pathogenic, pathogenic (dominant)',
#        'no classification for the single variant',
#        'likely pathogenic (dominant), pathogenic (dominant)',
#        'likely pathogenic (dominant), NA, VUS',
#        'likely pathogenic (dominant), NA, pathogenic',
#        'benign, likely benign, NA',
#        'likely pathogenic (dominant), NA, pathogenic, pathogenic (dominant), VUS',
#        'likely pathogenic (dominant), pathogenic, VUS',
#        'likely pathogenic, likely pathogenic (dominant), pathogenic',
#        'likely pathogenic, NA, pathogenic',
#        'likely pathogenic, likely pathogenic (dominant), NA, pathogenic, pathogenic (dominant)',
#        'likely pathogenic, NA, pathogenic, pathogenic (dominant), VUS',
#        'benign, likely benign, NA, pathogenic (dominant), VUS',
#        'likely benign, likely pathogenic', 'NA, pathogenic (dominant)',
#        'likely pathogenic, pathogenic, pathogenic (dominant), VUS',
#        'likely benign (!), NA, VUS', 'benign, pathogenic', 'VUS (!)',
#        'benign, likely benign, likely pathogenic, VUS',
#        'likely pathogenic (dominant), pathogenic, pathogenic (dominant), VUS',
#        'likely benign, likely pathogenic, likely pathogenic (dominant), VUS',
#        'likely benign, NA, pathogenic, VUS',
#        'likely benign (dominant), NA, VUS',
#        'likely benign, pathogenic, VUS',
#        'benign, likely benign, NA, VUS, VUS (!)',
#        'likely benign, likely pathogenic, VUS',
#        'benign, likely pathogenic (dominant), NA, VUS',
#        'benign, likely benign, NA, pathogenic, VUS',
#        'likely benign (dominant), VUS',
#        'benign, likely benign, likely pathogenic, NA, VUS',
#        'likely benign, pathogenic, pathogenic (dominant)',
#        'pathogenic (!)', 'benign, likely benign, pathogenic, VUS',
#        'likely benign, likely pathogenic, pathogenic, VUS',
#        'likely pathogenic (dominant), NA, pathogenic, VUS',
#        'likely pathogenic, NA, VUS',
#        'likely pathogenic (dominant), NA, pathogenic (dominant)',
#        'likely pathogenic (dominant), VUS',
#        'likely pathogenic, pathogenic (recessive)',
#        'likely pathogenic, likely pathogenic (dominant), NA',
#        'benign, likely benign, likely pathogenic',
#        'benign, likely pathogenic, NA, pathogenic, VUS',
#        'likely pathogenic, likely pathogenic (dominant), NA, VUS',
#        'benign, likely pathogenic (dominant), NA, pathogenic, pathogenic (dominant)',
#        'benign, likely benign, likely pathogenic, NA',
#        'likely benign, pathogenic', 'NA, pathogenic (!)',
#        'benign, likely benign, NA, pathogenic',
#        'NA, pathogenic (dominant), VUS',
#        'benign, likely pathogenic (dominant), NA, pathogenic, pathogenic (dominant), VUS',
#        'benign, pathogenic, pathogenic (dominant)',
#        'benign, likely pathogenic, NA, VUS',
#        'likely pathogenic, likely pathogenic (dominant), NA, pathogenic, VUS',
#        'likely pathogenic, likely pathogenic (dominant)',
#        'likely pathogenic, likely pathogenic (dominant), NA, pathogenic',
#        'likely pathogenic, likely pathogenic (dominant), NA, pathogenic (dominant), VUS',
#        'likely pathogenic (!), NA, VUS',
#        'likely benign, likely pathogenic, NA',
#        'likely benign, likely pathogenic, NA, VUS',
#        'pathogenic (dominant), pathogenic (recessive)',
#        'pathogenic, pathogenic (recessive)', 'pathogenic (recessive)',
#        'benign, likely benign, pathogenic (recessive), VUS',
#        'likely pathogenic, pathogenic, pathogenic (recessive), VUS',
#        'likely pathogenic (recessive)',
#        'benign, likely benign, pathogenic (recessive)',
#        'likely pathogenic (recessive), VUS',
#        'pathogenic (recessive), VUS', 'benign, pathogenic (recessive)',
#        'likely pathogenic, pathogenic, pathogenic (recessive)',
#        'likely pathogenic (recessive), pathogenic, pathogenic (recessive), VUS',
#        'likely pathogenic (recessive), pathogenic (recessive)',
#        'likely benign, pathogenic (recessive), VUS',
#        'pathogenic, pathogenic (recessive), VUS',
#        'likely benign, pathogenic (recessive)',
#        'likely pathogenic (recessive), pathogenic',
#        'benign, likely benign, pathogenic', 'NA, pathogenic (recessive)',
#        'likely pathogenic (recessive), NA, pathogenic, pathogenic (recessive), VUS',
#        'likely pathogenic (recessive), pathogenic, pathogenic (recessive)',
#        'likely pathogenic, NA, pathogenic (recessive)',
#        'NA, pathogenic, pathogenic (recessive)', 'likely pathogenic (!)',
#        'pathogenic, pathogenic (!), pathogenic (dominant)',
#        'Likely pathogenic, low penetrance',
#        'NA, pathogenic (recessive), VUS',
#        'NA, pathogenic, pathogenic (recessive), VUS'
# ]

# # Dictionary to store user choices
# label_mapping = {}

# # Mapping from user input to final label
# choice_labels = {
#     '1': 'pathogenic',
#     '2': 'benign',
#     '3': 'VUS'
# }

# for conseq in conseq_list:
#     while True:
#         print("\nCONSEQ value:")
#         print(f"  {conseq}")
#         print("Choose a label:")
#         print("  1) pathogenic")
#         print("  2) benign")
#         print("  3) VUS")
#         user_input = input("Enter 1, 2, or 3: ").strip()
#         if user_input in choice_labels:
#             label_mapping[conseq] = choice_labels[user_input]
#             break
#         else:
#             print("Invalid input. Please enter 1, 2, or 3.")

# # Once finished, display the results
# print("\nFinal CONSEQ labeling:")
# for conseq, label in label_mapping.items():
#     print(f"{conseq!r}: {label}")


import re

# Multi-line string of the final CONSEQ labeling
mapping_str = """
'Pathogenic': pathogenic
'benign': benign
'Benign': benign
'Uncertain significance': VUS
'benign, VUS': benign
'Likely benign': benign
'Conflicting classifications of pathogenicity': VUS
'VUS': VUS
'Benign/Likely benign': benign
'likely benign': benign
'benign, likely benign, VUS': benign
'-': VUS
'nan': VUS
'likely benign, VUS': benign
'pathogenic, VUS': pathogenic
'Likely pathogenic': pathogenic
'benign, NA': benign
'pathogenic (dominant)': pathogenic
'NA, VUS': VUS
'pathogenic, pathogenic (dominant)': pathogenic
'not provided': VUS
'NA, pathogenic, pathogenic (dominant)': pathogenic
'NA, pathogenic, VUS': pathogenic
'Pathogenic/Likely pathogenic': pathogenic
'pathogenic': pathogenic
'pathogenic, pathogenic (dominant), VUS': pathogenic
'likely benign, NA': benign
'likely pathogenic': pathogenic
'likely pathogenic (dominant), NA': pathogenic
'likely pathogenic, VUS': pathogenic
'likely pathogenic, NA': pathogenic
'likely pathogenic, NA, pathogenic, VUS': pathogenic
'likely pathogenic (dominant), pathogenic, pathogenic (dominant)': pathogenic
'benign, NA, VUS': benign
'likely pathogenic (dominant)': pathogenic
'NA, pathogenic': pathogenic
'likely pathogenic, NA, pathogenic, pathogenic (dominant)': pathogenic
'pathogenic (dominant), VUS': pathogenic
'likely pathogenic, likely pathogenic (dominant), NA, pathogenic, pathogenic (dominant), VUS': pathogenic
'likely benign, likely pathogenic (dominant), pathogenic, VUS': pathogenic
'benign, likely benign': benign
'benign, likely benign, NA, VUS': benign
'likely pathogenic, pathogenic, pathogenic (dominant)': pathogenic
'likely pathogenic (dominant), pathogenic': pathogenic
'likely benign, NA, VUS': benign
'likely pathogenic, pathogenic (dominant)': pathogenic
'likely pathogenic, pathogenic': pathogenic
'benign, NA, pathogenic, pathogenic (dominant), VUS': VUS
'NA, pathogenic, pathogenic (dominant), VUS': pathogenic
'benign, likely pathogenic': VUS
'likely pathogenic, pathogenic, VUS': pathogenic
'likely pathogenic (dominant), NA, pathogenic, pathogenic (dominant)': pathogenic
'no classification for the single variant': VUS
'likely pathogenic (dominant), pathogenic (dominant)': pathogenic
'likely pathogenic (dominant), NA, VUS': pathogenic
'likely pathogenic (dominant), NA, pathogenic': pathogenic
'benign, likely benign, NA': benign
'likely pathogenic (dominant), NA, pathogenic, pathogenic (dominant), VUS': pathogenic
'likely pathogenic (dominant), pathogenic, VUS': pathogenic
'likely pathogenic, likely pathogenic (dominant), pathogenic': pathogenic
'likely pathogenic, NA, pathogenic': pathogenic
'likely pathogenic, likely pathogenic (dominant), NA, pathogenic, pathogenic (dominant)': pathogenic
'likely pathogenic, NA, pathogenic, pathogenic (dominant), VUS': pathogenic
'benign, likely benign, NA, pathogenic (dominant), VUS': VUS
'likely benign, likely pathogenic': VUS
'NA, pathogenic (dominant)': pathogenic
'likely pathogenic, pathogenic, pathogenic (dominant), VUS': pathogenic
'likely benign (!), NA, VUS': benign
'benign, pathogenic': VUS
'VUS (!)': VUS
'benign, likely benign, likely pathogenic, VUS': VUS
'likely pathogenic (dominant), pathogenic, pathogenic (dominant), VUS': pathogenic
'likely benign, likely pathogenic, likely pathogenic (dominant), VUS': VUS
'likely benign, NA, pathogenic, VUS': VUS
'likely benign (dominant), NA, VUS': benign
'likely benign, pathogenic, VUS': VUS
'benign, likely benign, NA, VUS, VUS (!)': benign
'likely benign, likely pathogenic, VUS': VUS
'benign, likely pathogenic (dominant), NA, VUS': VUS
'benign, likely benign, NA, pathogenic, VUS': VUS
'likely benign (dominant), VUS': benign
'benign, likely benign, likely pathogenic, NA, VUS': VUS
'likely benign, pathogenic, pathogenic (dominant)': VUS
'pathogenic (!)': pathogenic
'benign, likely benign, pathogenic, VUS': VUS
'likely benign, likely pathogenic, pathogenic, VUS': VUS
'likely pathogenic (dominant), NA, pathogenic, VUS': pathogenic
'likely pathogenic, NA, VUS': pathogenic
'likely pathogenic (dominant), NA, pathogenic (dominant)': pathogenic
'likely pathogenic (dominant), VUS': pathogenic
'likely pathogenic, pathogenic (recessive)': pathogenic
'likely pathogenic, likely pathogenic (dominant), NA': pathogenic
'benign, likely benign, likely pathogenic': VUS
'benign, likely pathogenic, NA, pathogenic, VUS': VUS
'likely pathogenic, likely pathogenic (dominant), NA, VUS': pathogenic
'benign, likely pathogenic (dominant), NA, pathogenic, pathogenic (dominant)': VUS
'benign, likely benign, likely pathogenic, NA': VUS
'likely benign, pathogenic': VUS
'NA, pathogenic (!)': pathogenic
'benign, likely benign, NA, pathogenic': VUS
'NA, pathogenic (dominant), VUS': pathogenic
'benign, likely pathogenic (dominant), NA, pathogenic, pathogenic (dominant), VUS': VUS
'benign, pathogenic, pathogenic (dominant)': VUS
'benign, likely pathogenic, NA, VUS': VUS
'likely pathogenic, likely pathogenic (dominant), NA, pathogenic, VUS': pathogenic
'likely pathogenic, likely pathogenic (dominant)': pathogenic
'likely pathogenic, likely pathogenic (dominant), NA, pathogenic': pathogenic
'likely pathogenic, likely pathogenic (dominant), NA, pathogenic (dominant), VUS': pathogenic
'likely pathogenic (!), NA, VUS': pathogenic
'likely benign, likely pathogenic, NA': VUS
'likely benign, likely pathogenic, NA, VUS': VUS
'pathogenic (dominant), pathogenic (recessive)': pathogenic
'pathogenic, pathogenic (recessive)': pathogenic
'pathogenic (recessive)': pathogenic
'benign, likely benign, pathogenic (recessive), VUS': VUS
'likely pathogenic, pathogenic, pathogenic (recessive), VUS': pathogenic
'likely pathogenic (recessive)': pathogenic
'benign, likely benign, pathogenic (recessive)': VUS
'likely pathogenic (recessive), VUS': pathogenic
'pathogenic (recessive), VUS': pathogenic
'benign, pathogenic (recessive)': VUS
'likely pathogenic, pathogenic, pathogenic (recessive)': pathogenic
'likely pathogenic (recessive), pathogenic, pathogenic (recessive), VUS': pathogenic
'likely pathogenic (recessive), pathogenic (recessive)': pathogenic
'likely benign, pathogenic (recessive), VUS': VUS
'pathogenic, pathogenic (recessive), VUS': pathogenic
'likely benign, pathogenic (recessive)': VUS
'likely pathogenic (recessive), pathogenic': pathogenic
'benign, likely benign, pathogenic': VUS
'NA, pathogenic (recessive)': pathogenic
'likely pathogenic (recessive), NA, pathogenic, pathogenic (recessive), VUS': pathogenic
'likely pathogenic (recessive), pathogenic, pathogenic (recessive)': pathogenic
'likely pathogenic, NA, pathogenic (recessive)': pathogenic
'NA, pathogenic, pathogenic (recessive)': pathogenic
'likely pathogenic (!)': pathogenic
'pathogenic, pathogenic (!), pathogenic (dominant)': pathogenic
'Likely pathogenic, low penetrance': pathogenic
'NA, pathogenic (recessive), VUS': pathogenic
'NA, pathogenic, pathogenic (recessive), VUS': pathogenic
"""

# Parse into a dictionary
mapping = dict(re.findall(r"'([^']+)':\s*(\w+)", mapping_str))

# Display the resulting dictionary
print(mapping)

