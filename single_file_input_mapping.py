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

    `unmapped_data` is some data structure sent by
    the frontend.  The structure is known to the 
    developer since they specified the input element responsible
    for creating the data.  For example, a file chooser will send
    a list/array of primary keys.

    `id_list` is a list of WDL input "names"/ids that we are mapping
    to.  Note that the ordering is important.  Make sure the logic below
    matches the order in gui.json 

    '''
    unmapped_data = all_data[data_name]
    resource_pk = unmapped_data
    r = Resource.objects.get(pk=resource_pk)
    if (r.owner == user) or (user.is_staff):
        return {id_list[0]:r.path}
