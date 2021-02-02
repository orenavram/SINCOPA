def prepare(wd_path):

    import os
    import sys
    import shutil
    sys.path.insert(0, '/bioseq/pasa/')
    from CONSTANTS import RELOAD_TAGS, RESULT_WEBPAGE_NAME
    error_file = os.path.join(wd_path, 'error.txt')
    try:
        shutil.rmtree(f'{wd_path}/outputs')
        print(f'outputs were deleted!')
    except:
        print(f'No outputs dir to delete was found.')
    try:
        shutil.rmtree(f'{wd_path}/tmp')
        print(f'tmp files were deleted!')
    except:
        print(f'No tmp dir to delete was found.')

    html_path = os.path.join(wd_path, RESULT_WEBPAGE_NAME)
    with open(html_path) as f:
        html_content = ''
        for line in f:
            html_content += line
            if line.startswith('<!--result-->'):
                break

    html_content = html_content.replace('FAILED', 'RUNNING')
    html_content = html_content.replace('QUEUED', 'RUNNING')
    html_content = html_content.replace('FINISHED', 'RUNNING')
    if 'progress-bar-striped active' not in html_content:
        html_content = html_content.replace('progress-bar-striped', 'progress-bar-striped active')
    while f'<!--{RELOAD_TAGS}-->' in html_content:
        # maybe there is a double comment etc.. <!--<!--RELOAD_TAGS-->-->
        html_content = html_content.replace(f'<!--{RELOAD_TAGS}-->', RELOAD_TAGS)

    with open(html_path, 'w') as f:
        f.write(html_content)

    print(f'{html_path} was reverted!')


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('wd_path', help="A path to a working dir of pasa's job.\n"
                                        "E.g., /bioseq/data/results/sincopa/158875031326946667844691750504")
    args = parser.parse_args()

    prepare(args.wd_path)