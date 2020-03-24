from base.models import Resource
import os

def map_inputs(user, all_data, data_name, id_list):
    '''
    `user` is a User instance (or subclass).  This gives us
    the option of applying user-specific logic to the mapping.
    Since the code that calls this function does NOT know
    the structure of the input data, it cannot impose any logic
    such as filtering Resource objects for a particular user.
    Therefore we have to keep that information here

    `all_data` is some data structure sent by
    the frontend.  The structure is known to the 
    developer since they specified the input element responsible
    for creating the data.  For example, a file chooser will send
    a list/array of primary keys.

    `data_name` is how the data specific to this UI element is addressed
    in the `all_data` structure.

    `id_list` is a list of WDL input "names"/ids that we are mapping
    to.  Note that the ordering is important.  Make sure the logic below
    matches the order in gui.json 

    '''
    unmapped_data = all_data[data_name]

    r1_suffix = '_R1.fastq.gz'
    r2_suffix = '_R2.fastq.gz'
    r1_path_list = []
    r2_path_list = []
    for pk in unmapped_data:
        r = Resource.objects.get(pk=pk)
        if (r.owner == user) or (user.is_staff):
            if r.path.endswith(r1_suffix):
                r1_path_list.append(r.path)
            elif r.path.endswith(r2_suffix):
                r2_path_list.append(r.path)
            else:
                print('Skipping %s' % r.path)
        else:
            raise Exception('The user %s is not the owner of Resource with primary key %s.' % (user, pk))
    
    # now we have a list of files that had the correct naming scheme.
    # Need to check for pairing:
    r1_samples = [os.path.basename(x)[:-len(r1_suffix)] for x in r1_path_list]
    r2_samples = [os.path.basename(x)[:-len(r2_suffix)] for x in r2_path_list]
    r1_dict = dict(zip(r1_samples, r1_path_list))
    r2_dict = dict(zip(r2_samples, r2_path_list))

    sample_intersection = set(r1_samples).intersection(r2_samples)

    # now have the samples that have both R1 and R2.  Create the final map
    final_r1_list = []
    final_r2_list = []
    for s in sample_intersection:
        final_r1_list.append(r1_dict[s])
        final_r2_list.append(r2_dict[s])
    return {id_list[0]:final_r1_list, id_list[1]:final_r2_list}
