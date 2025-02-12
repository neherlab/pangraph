// Add index.html to request URLs without a file name in a CloudFront Functions viewer request event
// https://docs.aws.amazon.com/AmazonCloudFront/latest/DeveloperGuide/example_cloudfront_functions_url_rewrite_single_page_apps_section.html
async function handler(event) {
  const request = event.request;
  let uri = request.uri;

  // Check whether the URI ends with a slash (indicating it's likely a directory)
  if (uri.endsWith('/')) {
    request.uri += 'index.html';
  }

  // Check whether the URI is missing a file extension (indicating it might be a directory)
  else if (!uri.includes('.')) {
    request.uri += '/index.html';
  }

  return request;
}

if (typeof module !== 'undefined') {
  module.exports = handler;
}
