import os
import shutil

base_directory = 'mutatedRNAs-0.25'

file_to_find = 'pairagon_seeded.gtf'



def directory_has_file(directory, file_name):
    for root, dirs, files in os.walk(directory):
        if file_name in files:
            return True
    return False

subdirectories = [os.path.join(base_directory, d) for d in os.listdir(base_directory) if os.path.isdir(os.path.join(base_directory, d))]

kept = []
deleted = []


# Check each directory and remove it if it does not contain the specified file
for directory in subdirectories:

    # has_pairagon = directory_has_file(directory, file_to_find)
    # has_gmap = directory_has_file(directory, 'gmap.gff')

    # if not (has_pairagon and has_gmap):
    #     shutil.rmtree(directory)
    #     deleted.append(directory)
    #     # print(f"Deleted {directory}")
    # else:
    #     kept.append(directory)
    #     # print(f"Kept {directory}")


    file_path = os.path.join(directory, 'ans.txt')
    with open(file_path, 'r') as file:
        content = file.read()
    
    # Modify content if last character is '-'
    if content.endswith('-'):
        content_ = content[:-1] + '+'
    else:
        print(f"{directory} does not need changed")
        continue
    
    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.write(content_)




# kept.sort()
# deleted.sort()

# print("kept below")
# print("Completed 100 files")
# for folder in kept: 
#     print(folder)
# print("kept above^^")

# # print("deleted below")
# # for folder in deleted:
# #     print(folder)
# # print("deleted above^^")

# print(f"Kept: {len(kept)}")
# print(f"Deleted: {len(deleted)}")
