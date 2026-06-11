use crate::annotation::feature::Feature;
use crate::make_error;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_path::PathId;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};

/// Group annotation [`Feature`]s by the pangraph path they belong to.
///
/// Each feature's `seqid` — optionally remapped through `seqid_map` (the future `--seqid-map`
/// override) — is matched against pangraph path names. Every seqid that matches no path is
/// collected and reported together as a single hard error, rather than failing on the first
/// mismatch, since seqid/path-name mismatches are the most common failure mode.
pub fn match_features_to_paths(
  features: Vec<Feature>,
  graph: &Pangraph,
  seqid_map: &BTreeMap<String, String>,
) -> Result<BTreeMap<PathId, Vec<Feature>>, Report> {
  // Build a path-name -> PathId lookup once (paths without a name are skipped).
  let name_to_id: BTreeMap<&str, PathId> = graph
    .paths
    .values()
    .filter_map(|path| path.name().as_deref().map(|name| (name, path.id())))
    .collect();

  let mut grouped: BTreeMap<PathId, Vec<Feature>> = BTreeMap::new();
  let mut unmatched: BTreeSet<String> = BTreeSet::new();

  for feature in features {
    let path_name = seqid_map
      .get(&feature.seqid)
      .map_or(feature.seqid.as_str(), String::as_str);
    if let Some(&path_id) = name_to_id.get(path_name) {
      grouped.entry(path_id).or_default().push(feature);
    } else {
      unmatched.insert(feature.seqid.clone());
    }
  }

  if !unmatched.is_empty() {
    let names = unmatched.iter().cloned().collect::<Vec<_>>().join(", ");
    return make_error!(
      "Could not match {} annotation sequence id(s) to any pangraph path: {names}. \
       Check that annotation seqids correspond to the FASTA record names used to build the \
       graph, or provide an explicit seqid-to-path mapping.",
      unmatched.len()
    );
  }

  Ok(grouped)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pangraph::strand::Strand;
  use crate::utils::interval::Interval;
  use pretty_assertions::assert_eq;

  fn feat(seqid: &str) -> Feature {
    Feature {
      seqid: seqid.to_owned(),
      source: None,
      feature_type: "gene".to_owned(),
      interval: Interval::new(0, 10),
      strand: Some(Strand::Forward),
      id: None,
      name: None,
      attributes: vec![],
    }
  }

  fn test_graph() -> Result<Pangraph, Report> {
    Pangraph::from_path(&Some("../../data/test_graph.json"))
  }

  #[test]
  fn test_match_groups_features_by_path() -> Result<(), Report> {
    let graph = test_graph()?;
    let features = vec![feat("pKPC_CAV1321-45"), feat("pCAV1344-40"), feat("pKPC_CAV1321-45")];
    let grouped = match_features_to_paths(features, &graph, &BTreeMap::new())?;
    assert_eq!(grouped.len(), 2);
    let pid = graph.path_id_by_name("pKPC_CAV1321-45")?;
    assert_eq!(grouped[&pid].len(), 2);
    Ok(())
  }

  #[test]
  fn test_match_fails_on_unknown_seqid() -> Result<(), Report> {
    let graph = test_graph()?;
    let features = vec![feat("pCAV1344-40"), feat("does-not-exist"), feat("also-missing")];
    let err = match_features_to_paths(features, &graph, &BTreeMap::new()).unwrap_err();
    let msg = err.to_string();
    // Both unmatched seqids are reported together.
    assert!(msg.contains("does-not-exist"), "message was: {msg}");
    assert!(msg.contains("also-missing"), "message was: {msg}");
    Ok(())
  }

  #[test]
  fn test_match_uses_seqid_map_override() -> Result<(), Report> {
    let graph = test_graph()?;
    let mut seqid_map = BTreeMap::new();
    seqid_map.insert("weird_id".to_owned(), "pCAV1669-34".to_owned());
    let grouped = match_features_to_paths(vec![feat("weird_id")], &graph, &seqid_map)?;
    let pid = graph.path_id_by_name("pCAV1669-34")?;
    assert!(grouped.contains_key(&pid));
    Ok(())
  }
}
