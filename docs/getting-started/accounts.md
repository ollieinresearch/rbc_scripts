# Accounts and access

Everything runs on the clusters of the [Digital Research Alliance of
Canada](https://alliancecan.ca) (formerly Compute Canada). The steps below are a
checklist; the Alliance documentation linked from each step is the authority on the
details, which change from time to time.

1. **Alliance account.** Register in the [CCDB](https://ccdb.alliancecan.ca) as a student,
   with David Goluskin as your sponsor ([Apply for a CCDB
   account](https://docs.alliancecan.ca/wiki/Apply_for_a_CCDB_account)). David approves the
   request.
2. **Multifactor authentication.** Enrol a second factor; logins to the clusters require
   it ([Multifactor authentication](https://docs.alliancecan.ca/wiki/Multifactor_authentication)).
3. **SSH key.** Create a key on your laptop and upload the public key in the CCDB
   ([SSH keys](https://docs.alliancecan.ca/wiki/SSH_Keys)), then log in with
   `ssh <user>@<login address>`, the address being on each cluster's Alliance page
   ([Fir](https://docs.alliancecan.ca/wiki/Fir), [Trillium](https://docs.alliancecan.ca/wiki/Trillium)).
4. **Allocations.** Being sponsored by David gives you the group's default account
   `def-goluskin` (storage under `/project/def-goluskin/<user>`) and, while the group
   holds one, the competitive allocation `rrg-goluskin`, which the job scripts charge
   (`#SBATCH --account=rrg-goluskin`). Ask David which account new runs should charge.
5. **Trillium.** Trillium is run by SciNet; see the
   [Trillium page](https://docs.alliancecan.ca/wiki/Trillium) for whether access has to be
   requested separately in the CCDB.
6. **GitHub.** An account, so you can clone the repository and push changes.

```{admonition} Where things live
:class: note

Your files are in three places on every cluster: {term}`home` (small, backed up),
{term}`scratch` (large and fast, not backed up) and {term}`project space` (group storage,
backed up). Simulations run on scratch and are mirrored to home and project space
automatically; see {doc}`../user-guide/storage`.
```

Next: {doc}`install`.
